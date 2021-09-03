# Copyright 2021 by Michiel de Hoon.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.Align.fasta module."""
import unittest
import os

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Align.fasta_m8 import AlignmentIterator

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.emboss."
    ) from None


class TestFastaProtein(unittest.TestCase):

    query = Seq("MPMILGYWNVRGLTHPIRMLLEYTDSSYDEKRYTMGDAPDFDRSQWLNEKFKLGLDFPNLPYLIDGSHKITQSNAILRYLARKHHLDGETEEERIRADIVENQVMDTRMQLIMLCYNPDFEKQKPEFLKTIPEKMKLYSEFLGKRPWFAGDKVTYVDFLAYDILDQYRMFEPKCLDAFPNLRDFLARFEGLKKISAYMKSSRYIATPIFSKMAHWSNK")

    filename = os.path.join("Fasta", "protein_lib.fa")
    records = SeqIO.parse(filename, 'fasta')
    targets = {record.id: record.seq.upper() for record in records}

    def test_m8CB(self):
        # Alignment file obtained by running
        # fasta36 -q -m 8CB seq/mgstm1.aa seq/prot_test.lseg
        # in the fasta36 source distribution
        path = "Fasta/protein_m8CB.txt"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(
                alignments.commandline,
                "fasta36 -q -m 8CB seq/mgstm1.aa seq/prot_test.lseg",
            )
            alignments = list(alignments)
        self.assertEqual(len(alignments), 12)
        # sp|P10649|GSTM1_MOUSE   sp|P09488|GSTM1_HUMAN
        alignment = alignments[0]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 218)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 218)
        self.assertEqual(alignment.shape, (2, 218 + 0))
        self.assertEqual(alignment.sequences[0].id, "sp|P09488|GSTM1_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 48)
        self.assertAlmostEqual(alignment.annotations["evalue"], 6.1e-78)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 275.6)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[0, 218], [0, 218]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), len(target))
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0],
            "MPMILGYWDIRGLAHAIRLLLEYTDSSYEEKKYTMGDAPDYDRSQWLNEKFKLGLDFPNLPYLIDGAHKITQSNAILCYIARKHNLCGETEEEKIRVDILENQTMDNHMQLGMICYNPEFEKLKPKYLEELPEKLKLYSEFLGKRPWFAGNKITFVDFLVYDVLDLHRIFEPKCLDAFPNLKDFISRFEGLEKISAYMKSSRFLPRPVFSKMAVWGNK",
        )
        self.assertEqual(
            alignment[1],
            "MPMILGYWNVRGLTHPIRMLLEYTDSSYDEKRYTMGDAPDFDRSQWLNEKFKLGLDFPNLPYLIDGSHKITQSNAILRYLARKHHLDGETEEERIRADIVENQVMDTRMQLIMLCYNPDFEKQKPEFLKTIPEKMKLYSEFLGKRPWFAGDKVTYVDFLAYDILDQYRMFEPKCLDAFPNLRDFLARFEGLKKISAYMKSSRYIATPIFSKMAHWSNK",
        )
        # sp|P10649|GSTM1_MOUSE   sp|P00502|GSTA1_RAT
        alignment = alignments[1]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 205)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 205)
        self.assertEqual(alignment.shape, (2, 205 + 18))
        self.assertEqual(alignment.sequences[0].id, "sp|P00502|GSTA1_RAT")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 144)
        self.assertAlmostEqual(alignment.annotations["evalue"], 3.1e-13)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 60.7)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [
                            5,
                            33,
                            33,
                            46,
                            48,
                            58,
                            59,
                            62,
                            62,
                            101,
                            102,
                            125,
                            127,
                            142,
                            144,
                            218,
                        ],
                        [
                            3,
                            31,
                            40,
                            53,
                            53,
                            63,
                            63,
                            66,
                            67,
                            106,
                            106,
                            129,
                            129,
                            144,
                            144,
                            218,
                        ],
                    ]
                ),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 218)
        self.assertEqual(len(alignment.sequences[1].seq), 218)
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0],
            "VLHYFNARGRMECIRWLLAAAGVEFDEK---------FIQSPEDLEKLKKDGNLMFDQVPMVEIDG-MKLAQTRAILNYIATKYDLYGKDMKERALIDMYTEGILDLTEMIMQLVICPPDQKEAKTALAKDRTKNRYLPAFEKVLKSHGQDYLVGNRLTRVDIHLLELLLYVEEFDASLLTSFPLLKAFKSRISSLPNVKKFLQPGSQRKLPVDAKQIEEARK",
        )
        self.assertEqual(
            alignment[1],
            "ILGYWNVRGLTHPIRMLLEYTDSSYDEKRYTMGDAPDFDRSQWLNEKFKL--GLDFPNLPYL-IDGSHKITQSNAILRYLARKHHLDGETEEERIRADIVENQVMD-TRMQLIMLCYNPDFEKQKPEFLK--TIPEKMKLYSEFLGK--RPWFAGDKVTYVDFLAYDILDQYRMFEPKCLDAFPNLRDFLARFEGLKKISAYMKSSRYIATPIFSKMAHWSNK",
        )
        # sp|P10649|GSTM1_MOUSE   sp|P69905|HBA_HUMAN
        alignment = alignments[2]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 37)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 37)
        self.assertEqual(alignment.shape, (2, 37 + 2))
        self.assertEqual(alignment.sequences[0].id, "sp|P69905|HBA_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 27)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.19)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 20.9)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[35, 48, 48, 58, 59, 73], [176, 189, 190, 200, 200, 214]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 73)
        self.assertEqual(len(alignment.sequences[1].seq), 218)
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "SFPTTKTYFPHFD-LSHGSAQVKGHGKKVADALTNAVAH")
        self.assertEqual(alignment[1], "AFPNLRDFLARFEGLKKISAYMKS-SRYIATPIFSKMAH")
        # sp|P10649|GSTM1_MOUSE   sp|P00517|KAPCA_BOVIN
        alignment = alignments[3]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 70)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 70)
        self.assertEqual(alignment.shape, (2, 70 + 2))
        self.assertEqual(alignment.sequences[0].id, "sp|P00517|KAPCA_BOVIN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 53)
        self.assertAlmostEqual(alignment.annotations["evalue"], 1.4)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 19.2)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[228, 274, 276, 300], [136, 182, 182, 206]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 300)
        self.assertEqual(len(alignment.sequences[1].seq), 218)
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0],
            "IYEMAAGYPPFFADQPIQIYEKIVSGKVRFPSHFSSDLKDLLRNLLQVDLTKRFGNLKNGVNDIKNHKWFAT",
        )
        self.assertEqual(
            alignment[1],
            "LYSEFLGKRPWFAGDKVTYVDFLAYDILDQYRMFEPKCLDAFPNLR--DFLARFEGLKKISAYMKSSRYIAT",
        )
        # sp|P10649|GSTM1_MOUSE   sp|P14960|RBS_GUITH
        alignment = alignments[4]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 7)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 7)
        self.assertEqual(alignment.shape, (2, 7 + 0))
        self.assertEqual(alignment.sequences[0].id, "sp|P14960|RBS_GUITH")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 3)
        self.assertAlmostEqual(alignment.annotations["evalue"], 1.6)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 17.7)
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[46, 53], [6, 13]]))
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 53)
        self.assertEqual(len(alignment.sequences[1].seq), 218)
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "YWDLWGL")
        self.assertEqual(alignment[1], "YWNVRGL")
        # sp|P10649|GSTM1_MOUSE   sp|P01593|KV101_HUMAN
        alignment = alignments[5]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 40)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 40)
        self.assertEqual(alignment.shape, (2, 40 + 10))
        self.assertEqual(alignment.sequences[0].id, "sp|P01593|KV101_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 22)
        self.assertAlmostEqual(alignment.annotations["evalue"], 2.5)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 16.6)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [15, 29, 32, 40, 43, 47, 49, 51, 51, 58, 59, 64],
                        [149, 163, 163, 171, 171, 175, 175, 177, 178, 185, 185, 190],
                    ],
                ),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 64)
        self.assertEqual(len(alignment.sequences[1].seq), 218)
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0], "GDRVTITCQASQDINHYLNWYQQGPKKAPKILIYDA-SNLETGVPSRFSG"
        )
        self.assertEqual(
            alignment[1], "GDKVTYVDFLAYDI---LDQYRMFE---PKCL--DAFPNLRDFL-ARFEG"
        )
        # sp|P10649|GSTM1_MOUSE   sp|P99998|CYC_PANTR
        alignment = alignments[6]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 58)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 58)
        self.assertEqual(alignment.shape, (2, 58 + 10))
        self.assertEqual(alignment.sequences[0].id, "sp|P99998|CYC_PANTR")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 40)
        self.assertAlmostEqual(alignment.annotations["evalue"], 2.7)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 16.4)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [27, 47, 50, 58, 58, 68, 68, 73, 73, 82, 82, 88],
                        [128, 148, 148, 156, 157, 167, 168, 173, 175, 184, 187, 193],
                    ],
                ),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 88)
        self.assertEqual(len(alignment.sequences[1].seq), 218)
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0],
            "KTGPNLHGLFGRKTGQAPGYSYTAANKNKGI-IWGEDTLMEY-LENPK--KYIPGTKMI---FVGIKK",
        )
        self.assertEqual(
            alignment[1],
            "KTIPEKMKLYSEFLGKRPWF---AGDKVTYVDFLAYDILDQYRMFEPKCLDAFPNLRDFLARFEGLKK",
        )
        # sp|P10649|GSTM1_MOUSE   sp|P02585|TNNC2_HUMAN
        alignment = alignments[7]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 45)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 45)
        self.assertEqual(alignment.shape, (2, 45 + 9))
        self.assertEqual(alignment.sequences[0].id, "sp|P02585|TNNC2_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 31)
        self.assertAlmostEqual(alignment.annotations["evalue"], 3.9)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 16.4)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [[12, 26, 26, 36, 38, 49, 52, 62], [43, 57, 61, 71, 71, 82, 82, 92]]
                ),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 62)
        self.assertEqual(len(alignment.sequences[1].seq), 218)
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0], "SEEMIAEFKAAFDM----FDADGGGDISVKELGTVMRMLGQTPTKEELDAIIEE"
        )
        self.assertEqual(
            alignment[1], "SQWLNEKFKLGLDFPNLPYLIDGSHKIT--QSNAILRYLAR---KHHLDGETEE"
        )
        # sp|P10649|GSTM1_MOUSE   sp|P60615|NXL1A_BUNMU
        alignment = alignments[8]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 9)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 9)
        self.assertEqual(alignment.shape, (2, 9 + 2))
        self.assertEqual(alignment.sequences[0].id, "sp|P60615|NXL1A_BUNMU")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 3)
        self.assertAlmostEqual(alignment.annotations["evalue"], 4.2)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 15.6)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[85, 86, 86, 89, 89, 94], [114, 115, 116, 119, 120, 125]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 94)
        self.assertEqual(len(alignment.sequences[1].seq), 218)
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "C-NPH-PKQRP")
        self.assertEqual(alignment[1], "CYNPDFEKQKP")
        # sp|P10649|GSTM1_MOUSE   sp|P00193|FER_PEPAS
        alignment = alignments[9]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 4)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 4)
        self.assertEqual(alignment.shape, (2, 4 + 0))
        self.assertEqual(alignment.sequences[0].id, "sp|P00193|FER_PEPAS")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 2)
        self.assertAlmostEqual(alignment.annotations["evalue"], 4.2)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 14.7)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[14, 18], [170, 174]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 18)
        self.assertEqual(len(alignment.sequences[1].seq), 218)
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "KPEC")
        self.assertEqual(alignment[1], "EPKC")
        # sp|P10649|GSTM1_MOUSE   sp|P03435|HEMA_I75A3
        alignment = alignments[10]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 38)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 38)
        self.assertEqual(alignment.shape, (2, 38 + 1))
        self.assertEqual(alignment.sequences[0].id, "sp|P03435|HEMA_I75A3")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 27)
        self.assertAlmostEqual(alignment.annotations["evalue"], 4.7)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 17.9)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[398, 408, 409, 437], [73, 83, 83, 111]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 437)
        self.assertEqual(len(alignment.sequences[1].seq), 218)
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "NRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDL")
        self.assertEqual(alignment[1], "NAILRYLARK-HHLDGETEEERIRADIVENQVMDTRMQL")
        # sp|P10649|GSTM1_MOUSE   sp|P01834|IGKC_HUMAN
        alignment = alignments[11]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 67)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 67)
        self.assertEqual(alignment.shape, (2, 67 + 4))
        self.assertEqual(alignment.sequences[0].id, "sp|P01834|IGKC_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 54)
        self.assertAlmostEqual(alignment.annotations["evalue"], 6.3)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 14.9)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[11, 25, 27, 67, 69, 82], [57, 71, 71, 111, 111, 124]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 82)
        self.assertEqual(len(alignment.sequences[1].seq), 218)
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0],
            "PSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHK",
        )
        self.assertEqual(
            alignment[1],
            "PNLPYLIDGSHKIT--QSNAILRYLARKHHLDGETEEERIRADIVENQVMDTRMQL--IMLCYNPDFEKQK",
        )

    def test_m8CC(self):
        # Alignment file obtained by running
        # fasta36 -q -m 8CB seq/mgstm1.aa seq/prot_test.lseg
        # in the fasta36 source distribution
        path = "Fasta/protein_m8CC.txt"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(
                alignments.commandline,
                "fasta36 -q -m 8CC seq/mgstm1.aa seq/prot_test.lseg",
            )
            alignments = list(alignments)
        self.assertEqual(len(alignments), 12)
        # sp|P10649|GSTM1_MOUSE   sp|P09488|GSTM1_HUMAN
        alignment = alignments[0]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 218)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 218)
        self.assertEqual(alignment.shape, (2, 218 + 0))
        self.assertEqual(alignment.sequences[0].id, "sp|P09488|GSTM1_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 48)
        self.assertAlmostEqual(alignment.annotations["evalue"], 7.6e-83)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 291.9)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[0, 218], [0, 218]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), len(target))
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0],
            "MPMILGYWDIRGLAHAIRLLLEYTDSSYEEKKYTMGDAPDYDRSQWLNEKFKLGLDFPNLPYLIDGAHKITQSNAILCYIARKHNLCGETEEEKIRVDILENQTMDNHMQLGMICYNPEFEKLKPKYLEELPEKLKLYSEFLGKRPWFAGNKITFVDFLVYDVLDLHRIFEPKCLDAFPNLKDFISRFEGLEKISAYMKSSRFLPRPVFSKMAVWGNK",
        )
        self.assertEqual(
            alignment[1],
            "MPMILGYWNVRGLTHPIRMLLEYTDSSYDEKRYTMGDAPDFDRSQWLNEKFKLGLDFPNLPYLIDGSHKITQSNAILRYLARKHHLDGETEEERIRADIVENQVMDTRMQLIMLCYNPDFEKQKPEFLKTIPEKMKLYSEFLGKRPWFAGDKVTYVDFLAYDILDQYRMFEPKCLDAFPNLRDFLARFEGLKKISAYMKSSRYIATPIFSKMAHWSNK",
        )
        # sp|P10649|GSTM1_MOUSE   sp|P00502|GSTA1_RAT
        alignment = alignments[1]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 205)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 205)
        self.assertEqual(alignment.shape, (2, 205 + 18))
        self.assertEqual(alignment.sequences[0].id, "sp|P00502|GSTA1_RAT")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 144)
        self.assertAlmostEqual(alignment.annotations["evalue"], 4.4e-14)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 63.5)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [
                            5,
                            33,
                            33,
                            46,
                            48,
                            58,
                            59,
                            62,
                            62,
                            101,
                            102,
                            125,
                            127,
                            142,
                            144,
                            218,
                        ],
                        [
                            3,
                            31,
                            40,
                            53,
                            53,
                            63,
                            63,
                            66,
                            67,
                            106,
                            106,
                            129,
                            129,
                            144,
                            144,
                            218,
                        ],
                    ]
                ),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 218)
        self.assertEqual(len(alignment.sequences[1].seq), 218)
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0],
            "VLHYFNARGRMECIRWLLAAAGVEFDEK---------FIQSPEDLEKLKKDGNLMFDQVPMVEIDG-MKLAQTRAILNYIATKYDLYGKDMKERALIDMYTEGILDLTEMIMQLVICPPDQKEAKTALAKDRTKNRYLPAFEKVLKSHGQDYLVGNRLTRVDIHLLELLLYVEEFDASLLTSFPLLKAFKSRISSLPNVKKFLQPGSQRKLPVDAKQIEEARK",
        )
        self.assertEqual(
            alignment[1],
            "ILGYWNVRGLTHPIRMLLEYTDSSYDEKRYTMGDAPDFDRSQWLNEKFKL--GLDFPNLPYL-IDGSHKITQSNAILRYLARKHHLDGETEEERIRADIVENQVMD-TRMQLIMLCYNPDFEKQKPEFLK--TIPEKMKLYSEFLGK--RPWFAGDKVTYVDFLAYDILDQYRMFEPKCLDAFPNLRDFLARFEGLKKISAYMKSSRYIATPIFSKMAHWSNK",
        )
        # sp|P10649|GSTM1_MOUSE   sp|P69905|HBA_HUMAN
        alignment = alignments[2]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 37)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 37)
        self.assertEqual(alignment.shape, (2, 37 + 2))
        self.assertEqual(alignment.sequences[0].id, "sp|P69905|HBA_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 27)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.15)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 21.2)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[35, 48, 48, 58, 59, 73], [176, 189, 190, 200, 200, 214]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 73)
        self.assertEqual(len(alignment.sequences[1].seq), 218)
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "SFPTTKTYFPHFD-LSHGSAQVKGHGKKVADALTNAVAH")
        self.assertEqual(alignment[1], "AFPNLRDFLARFEGLKKISAYMKS-SRYIATPIFSKMAH")
        # sp|P10649|GSTM1_MOUSE   sp|P00517|KAPCA_BOVIN
        alignment = alignments[3]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 70)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 70)
        self.assertEqual(alignment.shape, (2, 70 + 2))
        self.assertEqual(alignment.sequences[0].id, "sp|P00517|KAPCA_BOVIN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 53)
        self.assertAlmostEqual(alignment.annotations["evalue"], 1.2)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 19.4)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[228, 274, 276, 300], [136, 182, 182, 206]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 300)
        self.assertEqual(len(alignment.sequences[1].seq), 218)
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0],
            "IYEMAAGYPPFFADQPIQIYEKIVSGKVRFPSHFSSDLKDLLRNLLQVDLTKRFGNLKNGVNDIKNHKWFAT",
        )
        self.assertEqual(
            alignment[1],
            "LYSEFLGKRPWFAGDKVTYVDFLAYDILDQYRMFEPKCLDAFPNLR--DFLARFEGLKKISAYMKSSRYIAT",
        )
        # sp|P10649|GSTM1_MOUSE   sp|P14960|RBS_GUITH
        alignment = alignments[4]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 7)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 7)
        self.assertEqual(alignment.shape, (2, 7 + 0))
        self.assertEqual(alignment.sequences[0].id, "sp|P14960|RBS_GUITH")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 3)
        self.assertAlmostEqual(alignment.annotations["evalue"], 1.5)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 17.8)
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[46, 53], [6, 13]]))
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 53)
        self.assertEqual(len(alignment.sequences[1].seq), 218)
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "YWDLWGL")
        self.assertEqual(alignment[1], "YWNVRGL")
        # sp|P10649|GSTM1_MOUSE   sp|P01593|KV101_HUMAN
        alignment = alignments[5]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 40)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 40)
        self.assertEqual(alignment.shape, (2, 40 + 10))
        self.assertEqual(alignment.sequences[0].id, "sp|P01593|KV101_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 22)
        self.assertAlmostEqual(alignment.annotations["evalue"], 2.4)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 16.7)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [15, 29, 32, 40, 43, 47, 49, 51, 51, 58, 59, 64],
                        [149, 163, 163, 171, 171, 175, 175, 177, 178, 185, 185, 190],
                    ],
                ),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 64)
        self.assertEqual(len(alignment.sequences[1].seq), 218)
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0], "GDRVTITCQASQDINHYLNWYQQGPKKAPKILIYDA-SNLETGVPSRFSG"
        )
        self.assertEqual(
            alignment[1], "GDKVTYVDFLAYDI---LDQYRMFE---PKCL--DAFPNLRDFL-ARFEG"
        )
        # sp|P10649|GSTM1_MOUSE   sp|P99998|CYC_PANTR
        alignment = alignments[6]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 58)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 58)
        self.assertEqual(alignment.shape, (2, 58 + 10))
        self.assertEqual(alignment.sequences[0].id, "sp|P99998|CYC_PANTR")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 40)
        self.assertAlmostEqual(alignment.annotations["evalue"], 2.7)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 16.5)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [27, 47, 50, 58, 58, 68, 68, 73, 73, 82, 82, 88],
                        [128, 148, 148, 156, 157, 167, 168, 173, 175, 184, 187, 193],
                    ],
                ),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 88)
        self.assertEqual(len(alignment.sequences[1].seq), 218)
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0],
            "KTGPNLHGLFGRKTGQAPGYSYTAANKNKGI-IWGEDTLMEY-LENPK--KYIPGTKMI---FVGIKK",
        )
        self.assertEqual(
            alignment[1],
            "KTIPEKMKLYSEFLGKRPWF---AGDKVTYVDFLAYDILDQYRMFEPKCLDAFPNLRDFLARFEGLKK",
        )
        # sp|P10649|GSTM1_MOUSE   sp|P02585|TNNC2_HUMAN
        alignment = alignments[7]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 45)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 45)
        self.assertEqual(alignment.shape, (2, 45 + 9))
        self.assertEqual(alignment.sequences[0].id, "sp|P02585|TNNC2_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 31)
        self.assertAlmostEqual(alignment.annotations["evalue"], 3.8)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 16.5)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [[12, 26, 26, 36, 38, 49, 52, 62], [43, 57, 61, 71, 71, 82, 82, 92]]
                ),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 62)
        self.assertEqual(len(alignment.sequences[1].seq), 218)
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0], "SEEMIAEFKAAFDM----FDADGGGDISVKELGTVMRMLGQTPTKEELDAIIEE"
        )
        self.assertEqual(
            alignment[1], "SQWLNEKFKLGLDFPNLPYLIDGSHKIT--QSNAILRYLAR---KHHLDGETEE"
        )
        # sp|P10649|GSTM1_MOUSE   sp|P60615|NXL1A_BUNMU
        alignment = alignments[8]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 9)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 9)
        self.assertEqual(alignment.shape, (2, 9 + 2))
        self.assertEqual(alignment.sequences[0].id, "sp|P60615|NXL1A_BUNMU")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 3)
        self.assertAlmostEqual(alignment.annotations["evalue"], 4.2)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 15.6)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[85, 86, 86, 89, 89, 94], [114, 115, 116, 119, 120, 125]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 94)
        self.assertEqual(len(alignment.sequences[1].seq), 218)
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "C-NPH-PKQRP")
        self.assertEqual(alignment[1], "CYNPDFEKQKP")
        # sp|P10649|GSTM1_MOUSE   sp|P03435|HEMA_I75A3
        alignment = alignments[9]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 38)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 38)
        self.assertEqual(alignment.shape, (2, 38 + 1))
        self.assertEqual(alignment.sequences[0].id, "sp|P03435|HEMA_I75A3")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 27)
        self.assertAlmostEqual(alignment.annotations["evalue"], 4.4)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 18.1)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[398, 408, 409, 437], [73, 83, 83, 111]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 437)
        self.assertEqual(len(alignment.sequences[1].seq), 218)
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "NRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDL")
        self.assertEqual(alignment[1], "NAILRYLARK-HHLDGETEEERIRADIVENQVMDTRMQL")
        # sp|P10649|GSTM1_MOUSE   sp|P00193|FER_PEPAS
        alignment = alignments[10]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 4)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 4)
        self.assertEqual(alignment.shape, (2, 4 + 0))
        self.assertEqual(alignment.sequences[0].id, "sp|P00193|FER_PEPAS")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 2)
        self.assertAlmostEqual(alignment.annotations["evalue"], 4.4)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 14.7)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[14, 18], [170, 174]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 18)
        self.assertEqual(len(alignment.sequences[1].seq), 218)
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "KPEC")
        self.assertEqual(alignment[1], "EPKC")
        # sp|P10649|GSTM1_MOUSE   sp|P01834|IGKC_HUMAN
        alignment = alignments[11]
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 67)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 67)
        self.assertEqual(alignment.shape, (2, 67 + 4))
        self.assertEqual(alignment.sequences[0].id, "sp|P01834|IGKC_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 54)
        self.assertAlmostEqual(alignment.annotations["evalue"], 6.4)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 14.9)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[11, 25, 27, 67, 69, 82], [57, 71, 71, 111, 111, 124]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 82)
        self.assertEqual(len(alignment.sequences[1].seq), 218)
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0],
            "PSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHK",
        )
        self.assertEqual(
            alignment[1],
            "PNLPYLIDGSHKIT--QSNAILRYLARKHHLDGETEEERIRADIVENQVMDTRMQL--IMLCYNPDFEKQK",
        )


class TestFastaNucleotide(unittest.TestCase):

    query = Seq("ATGCCTATGATACTGGGATACTGGAACGTCCGCGGACTGACACACCCGATCCGCATGCTCCTGGAATACACAGACTCAAGCTATGATGAGAAGAGATACACCATGGGTGACGCTCCCGACTTTGACAGAAGCCAGTGGCTGAATGAGAAGTTCAAGCTGGGCCTGGACTTTCCCAATCTGCCTTACTTGATCGATGGATCACACAAGATCACCCAGAGCAATGCCATCCTGCGCTACCTTGCCCGAAAGCACCACCTGGATGGAGAGACAGAGGAGGAGAGGATCCGTGCAGACATTGTGGAGAACCAGGTCATGGACACCCGCATGCAGCTCATCATGCTCTGTTACAACCCTGACTTTGAGAAGCAGAAGCCAGAGTTCTTGAAGACCATCCCTGAGAAAATGAAGCTCTACTCTGAGTTCCTGGGCAAGAGGCCATGGTTTGCAGGGGACAAGGTCACCTATGTGGATTTCCTTGCTTATGACATTCTTGACCAGTACCGTATGTTTGAGCCCAAGTGCCTGGACGCCTTCCCAAACCTGAGGGACTTCCTGGCCCGCTTCGAGGGCCTCAAGAAGATCTCTGCCTACATGAAGAGTAGCCGCTACATCGCAACACCTATATTTTCAAAGATGGCCCACTGGAGTAACAAGTAG")

    filename = os.path.join("Fasta", "nucleotide_lib.fa")
    records = SeqIO.parse(filename, 'fasta')
    targets = {record.id: record.seq.upper() for record in records}



if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
