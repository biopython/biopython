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

    query = Seq(
        "MPMILGYWNVRGLTHPIRMLLEYTDSSYDEKRYTMGDAPDFDRSQWLNEKFKLGLDFPNLPYLIDGSHKITQSNAILRYLARKHHLDGETEEERIRADIVENQVMDTRMQLIMLCYNPDFEKQKPEFLKTIPEKMKLYSEFLGKRPWFAGDKVTYVDFLAYDILDQYRMFEPKCLDAFPNLRDFLARFEGLKKISAYMKSSRYIATPIFSKMAHWSNK"
    )

    filename = os.path.join("Fasta", "protein_lib.fa")
    records = SeqIO.parse(filename, "fasta")
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
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[0].seq), 218)
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
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
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
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "SFPTTKTYFPHFD-LSHGSAQVKGHGKKVADALTNAVAH")
        self.assertEqual(alignment[1], "AFPNLRDFLARFEGLKKISAYMKS-SRYIATPIFSKMAH")
        # sp|P10649|GSTM1_MOUSE   sp|P00517|KAPCA_BOVIN
        alignment = alignments[3]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
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
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "YWDLWGL")
        self.assertEqual(alignment[1], "YWNVRGL")
        # sp|P10649|GSTM1_MOUSE   sp|P01593|KV101_HUMAN
        alignment = alignments[5]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
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
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
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
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
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
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "C-NPH-PKQRP")
        self.assertEqual(alignment[1], "CYNPDFEKQKP")
        # sp|P10649|GSTM1_MOUSE   sp|P00193|FER_PEPAS
        alignment = alignments[9]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "KPEC")
        self.assertEqual(alignment[1], "EPKC")
        # sp|P10649|GSTM1_MOUSE   sp|P03435|HEMA_I75A3
        alignment = alignments[10]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "NRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDL")
        self.assertEqual(alignment[1], "NAILRYLARK-HHLDGETEEERIRADIVENQVMDTRMQL")
        # sp|P10649|GSTM1_MOUSE   sp|P01834|IGKC_HUMAN
        alignment = alignments[11]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
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
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[0].seq), 218)
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
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
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
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "SFPTTKTYFPHFD-LSHGSAQVKGHGKKVADALTNAVAH")
        self.assertEqual(alignment[1], "AFPNLRDFLARFEGLKKISAYMKS-SRYIATPIFSKMAH")
        # sp|P10649|GSTM1_MOUSE   sp|P00517|KAPCA_BOVIN
        alignment = alignments[3]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
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
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "YWDLWGL")
        self.assertEqual(alignment[1], "YWNVRGL")
        # sp|P10649|GSTM1_MOUSE   sp|P01593|KV101_HUMAN
        alignment = alignments[5]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
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
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
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
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
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
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "C-NPH-PKQRP")
        self.assertEqual(alignment[1], "CYNPDFEKQKP")
        # sp|P10649|GSTM1_MOUSE   sp|P03435|HEMA_I75A3
        alignment = alignments[9]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "NRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDL")
        self.assertEqual(alignment[1], "NAILRYLARK-HHLDGETEEERIRADIVENQVMDTRMQL")
        # sp|P10649|GSTM1_MOUSE   sp|P00193|FER_PEPAS
        alignment = alignments[10]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "KPEC")
        self.assertEqual(alignment[1], "EPKC")
        # sp|P10649|GSTM1_MOUSE   sp|P01834|IGKC_HUMAN
        alignment = alignments[11]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/prot_test.lseg")
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
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
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

    query = Seq(
        "ATGCCTATGATACTGGGATACTGGAACGTCCGCGGACTGACACACCCGATCCGCATGCTCCTGGAATACACAGACTCAAGCTATGATGAGAAGAGATACACCATGGGTGACGCTCCCGACTTTGACAGAAGCCAGTGGCTGAATGAGAAGTTCAAGCTGGGCCTGGACTTTCCCAATCTGCCTTACTTGATCGATGGATCACACAAGATCACCCAGAGCAATGCCATCCTGCGCTACCTTGCCCGAAAGCACCACCTGGATGGAGAGACAGAGGAGGAGAGGATCCGTGCAGACATTGTGGAGAACCAGGTCATGGACACCCGCATGCAGCTCATCATGCTCTGTTACAACCCTGACTTTGAGAAGCAGAAGCCAGAGTTCTTGAAGACCATCCCTGAGAAAATGAAGCTCTACTCTGAGTTCCTGGGCAAGAGGCCATGGTTTGCAGGGGACAAGGTCACCTATGTGGATTTCCTTGCTTATGACATTCTTGACCAGTACCGTATGTTTGAGCCCAAGTGCCTGGACGCCTTCCCAAACCTGAGGGACTTCCTGGCCCGCTTCGAGGGCCTCAAGAAGATCTCTGCCTACATGAAGAGTAGCCGCTACATCGCAACACCTATATTTTCAAAGATGGCCCACTGGAGTAACAAGTAG"
    )

    filename = os.path.join("Fasta", "nucleotide_lib.fa")
    records = SeqIO.parse(filename, "fasta")
    targets = {record.id: record.seq.upper() for record in records}

    def test_m8CB(self):
        # Alignment file obtained by running
        # fasta36 -m 8CB seq/mgstm1.nt seq/gst.nlib
        # in the fasta36 source distribution
        path = "Fasta/nucleotide_m8CB.txt"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(
                alignments.commandline,
                "fasta36 -m 8CB seq/mgstm1.nt seq/gst.nlib",
            )
            alignments = list(alignments)
        self.assertEqual(len(alignments), 12)
        # pGT875   pGT875
        alignment = alignments[0]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 657)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 657)
        self.assertEqual(alignment.shape, (2, 657 + 0))
        self.assertEqual(alignment.sequences[0].id, "pGT875")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 0)
        self.assertAlmostEqual(alignment.annotations["evalue"], 3.6e-194)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 666.0)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[37, 694], [0, 657]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 694)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0],
            "ATGCCTATGATACTGGGATACTGGAACGTCCGCGGACTGACACACCCGATCCGCATGCTCCTGGAATACACAGACTCAAGCTATGATGAGAAGAGATACACCATGGGTGACGCTCCCGACTTTGACAGAAGCCAGTGGCTGAATGAGAAGTTCAAGCTGGGCCTGGACTTTCCCAATCTGCCTTACTTGATCGATGGATCACACAAGATCACCCAGAGCAATGCCATCCTGCGCTACCTTGCCCGAAAGCACCACCTGGATGGAGAGACAGAGGAGGAGAGGATCCGTGCAGACATTGTGGAGAACCAGGTCATGGACACCCGCATGCAGCTCATCATGCTCTGTTACAACCCTGACTTTGAGAAGCAGAAGCCAGAGTTCTTGAAGACCATCCCTGAGAAAATGAAGCTCTACTCTGAGTTCCTGGGCAAGAGGCCATGGTTTGCAGGGGACAAGGTCACCTATGTGGATTTCCTTGCTTATGACATTCTTGACCAGTACCGTATGTTTGAGCCCAAGTGCCTGGACGCCTTCCCAAACCTGAGGGACTTCCTGGCCCGCTTCGAGGGCCTCAAGAAGATCTCTGCCTACATGAAGAGTAGCCGCTACATCGCAACACCTATATTTTCAAAGATGGCCCACTGGAGTAACAAGTAG",
        )
        self.assertEqual(
            alignment[1],
            "ATGCCTATGATACTGGGATACTGGAACGTCCGCGGACTGACACACCCGATCCGCATGCTCCTGGAATACACAGACTCAAGCTATGATGAGAAGAGATACACCATGGGTGACGCTCCCGACTTTGACAGAAGCCAGTGGCTGAATGAGAAGTTCAAGCTGGGCCTGGACTTTCCCAATCTGCCTTACTTGATCGATGGATCACACAAGATCACCCAGAGCAATGCCATCCTGCGCTACCTTGCCCGAAAGCACCACCTGGATGGAGAGACAGAGGAGGAGAGGATCCGTGCAGACATTGTGGAGAACCAGGTCATGGACACCCGCATGCAGCTCATCATGCTCTGTTACAACCCTGACTTTGAGAAGCAGAAGCCAGAGTTCTTGAAGACCATCCCTGAGAAAATGAAGCTCTACTCTGAGTTCCTGGGCAAGAGGCCATGGTTTGCAGGGGACAAGGTCACCTATGTGGATTTCCTTGCTTATGACATTCTTGACCAGTACCGTATGTTTGAGCCCAAGTGCCTGGACGCCTTCCCAAACCTGAGGGACTTCCTGGCCCGCTTCGAGGGCCTCAAGAAGATCTCTGCCTACATGAAGAGTAGCCGCTACATCGCAACACCTATATTTTCAAAGATGGCCCACTGGAGTAACAAGTAG",
        )
        # pGT875   RABGLTR
        alignment = alignments[1]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 646)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 646)
        self.assertEqual(alignment.shape, (2, 646 + 0))
        self.assertEqual(alignment.sequences[0].id, "RABGLTR")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 135)
        self.assertAlmostEqual(alignment.annotations["evalue"], 1.9e-118)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 414.4)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[33, 679], [0, 646]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 679)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0],
            "ATGCCCATGACGCTGGGTTACTGGGACGTCCGTGGGCTGGCTCTGCCAATCCGCATGCTCCTGGAATACACGGACACCAGCTATGAGGAAAAGAAATACACCATGGGGGATGCTCCCAACTATGACCAAAGCAAGTGGCTGAGTGAGAAGTTCACCCTGGGCCTGGACTTTCCCAATCTGCCCTACCTAATTGATGGGACTCACAAGCTCACGCAGAGCAACGCCATCCTGCGCTACCTGGCCCGCAAGCACGGCCTGTGTGGGGAGACGGAAGAGGAGAGGATTCGCGTGGACATTCTGGAGAATCAGCTGATGGACAACCGCTTCCAACTTGTAAACGTCTGCTACAGTCCCGACTTTGAGAAGCTCAAGCCCGAGTACCTGAAGGGGCTCCCTGAGAAGCTGCAGCTGTACTCGCAGTTCCTGGGAAGCCTCCCCTGGTTCGCAGGGGACAAGATCACCTTCGCCGATTTCCTTGTCTACGACGTTCTTGACCAGAACCGGATATTTGTGCCTGGGTGCCTGGACGCGTTCCCAAACCTGAAGGACTTTCATGTCCGCTTTGAGGGCCTGCCGAAGATCTCTGCCTACATGAAGTCCAGCCGCTTTATCCGAGTCCCTGTGTTTTTAAAGAAGGCCACGTGGA",
        )
        self.assertEqual(
            alignment[1],
            "ATGCCTATGATACTGGGATACTGGAACGTCCGCGGACTGACACACCCGATCCGCATGCTCCTGGAATACACAGACTCAAGCTATGATGAGAAGAGATACACCATGGGTGACGCTCCCGACTTTGACAGAAGCCAGTGGCTGAATGAGAAGTTCAAGCTGGGCCTGGACTTTCCCAATCTGCCTTACTTGATCGATGGATCACACAAGATCACCCAGAGCAATGCCATCCTGCGCTACCTTGCCCGAAAGCACCACCTGGATGGAGAGACAGAGGAGGAGAGGATCCGTGCAGACATTGTGGAGAACCAGGTCATGGACACCCGCATGCAGCTCATCATGCTCTGTTACAACCCTGACTTTGAGAAGCAGAAGCCAGAGTTCTTGAAGACCATCCCTGAGAAAATGAAGCTCTACTCTGAGTTCCTGGGCAAGAGGCCATGGTTTGCAGGGGACAAGGTCACCTATGTGGATTTCCTTGCTTATGACATTCTTGACCAGTACCGTATGTTTGAGCCCAAGTGCCTGGACGCCTTCCCAAACCTGAGGGACTTCCTGGCCCGCTTCGAGGGCCTCAAGAAGATCTCTGCCTACATGAAGAGTAGCCGCTACATCGCAACACCTATATTTTCAAAGATGGCCCACTGGA",
        )
        # pGT875   BTGST
        alignment = alignments[2]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 413)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 413)
        self.assertEqual(alignment.shape, (2, 413 + 21))
        self.assertEqual(alignment.sequences[0].id, "BTGST")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 167)
        self.assertAlmostEqual(alignment.annotations["evalue"], 1.2e-07)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 46.4)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [
                            227,
                            376,
                            376,
                            383,
                            384,
                            401,
                            401,
                            461,
                            466,
                            472,
                            473,
                            486,
                            488,
                            501,
                            505,
                            535,
                            537,
                            543,
                            543,
                            655,
                        ],
                        [
                            175,
                            324,
                            325,
                            332,
                            332,
                            349,
                            352,
                            412,
                            412,
                            418,
                            418,
                            431,
                            431,
                            444,
                            444,
                            474,
                            474,
                            480,
                            482,
                            594,
                        ],
                    ]
                ),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 655)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0],
            "AGCTCCCCAAGTTCCAGGACGGAGACCTCACGCTGTACCAGTCCAATGCCATCCTGCGGCACCTGGGCCGCACCCTCGGGCTGTATGGGAAGGACCAGCAGGAGGCGGCCCTGGTGGACATGGTGAATGACGGTGTAGAGGACCTTCGC-TGCAAATACGTCTCCCTCATTTACA---CCAACTACGAGGCGGGCAAGGAGGACTATGTGAAGGCGCTGCCCCAGCACCTGAAGCCTTTCGAGACCCTGCTGTCCCAGAACAAGGGTGGCCAGGCCTTCATCGTGGGCGACCAGATCTCCTTTGCGGACTACAACCTGCT--GGACCTGCTTCGGATTCACCAGGTCCTGGCCCCCAGCTGTCTGGACTCCTTCCCCCTGCTCTCAGCCTACGTGGCCCGTCTCAACTCCCGGCCCAAGCTCAAGGCCTTCCTG",
        )
        self.assertEqual(
            alignment[1],
            "ATCTGCCTTACTTGATCGATGGATCACACAAGATCACCCAGAGCAATGCCATCCTGCGCTACCTTGCCCGAAAGCACCACCTGGATGGAGAGACAGAGGAGGAGAGGATCCGTGCAGACATTGTGGAGAACCAGGTCATGGACACCCGCATGCAGCT-CATCATGCTCTGTTACAACCCTGACTTTGAGAAGCAGAAGCCAGAGTTCTTGAAGACCATCCCTGAGAAAATGAAGCTCT-----ACTCTG-AGTTCCTGGGCAA--GAGGCCATGGTTT----GCAGGGGACAAGGTCACCTATGTGGATTTC--CTTGCTTATGACATTCTTGACCAGTACCGTATGTTTGAGCCCAAGTGCCTGGACGCCTTCCCAAACCTGAGGGACTTCCTGGCCCGCTTCGAGGGCCTCAAGAAGATCTCTGCCTACATG",
        )
        # pGT875   RABGSTB
        alignment = alignments[3]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 127)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 127)
        self.assertEqual(alignment.shape, (2, 127 + 8))
        self.assertEqual(alignment.sequences[0].id, "RABGSTB")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 42)
        self.assertAlmostEqual(alignment.annotations["evalue"], 2.1e-07)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 45.6)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [
                            156,
                            171,
                            173,
                            190,
                            190,
                            201,
                            202,
                            260,
                            261,
                            272,
                            272,
                            279,
                            279,
                            287,
                        ],
                        [
                            158,
                            173,
                            173,
                            190,
                            192,
                            203,
                            203,
                            261,
                            261,
                            272,
                            273,
                            280,
                            281,
                            289,
                        ],
                    ]
                ),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 287)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0],
            "GGGTATTGATGTTCCAGCAAGTGCCCATGGTTGA--GATTGATGGGATGAAGCTGGTGCAGACCAGAGCCATTTTCAACTACATTGCAGACAAGCACAACCTGTATGGGAAAGACATA-AAGGAGA-GAGCCCTG",
        )
        self.assertEqual(
            alignment[1],
            "GGGCCTGGACTTTCC--CAATCTGCCTTACTTGATCGATGGATCACA-CAAGATCACCCAGAGCAATGCCATCCTGCGCTACCTTGCCCGAAAGCACCACCTGGAT-GGAGAGACAGAGGAGGAGAGGATCCGTG",
        )
        # pGT875   OCDHPR
        alignment = alignments[4]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 23)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 23)
        self.assertEqual(alignment.shape, (2, 23 + 1))
        self.assertEqual(alignment.sequences[0].id, "OCDHPR")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 2)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.0092)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 30.1)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[2302, 2319, 2319, 2325], [265, 282, 283, 289]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 2325)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "AGACAGAGGAGGAGAAG-TCTGTG")
        self.assertEqual(alignment[1], "AGACAGAGGAGGAGAGGATCCGTG")
        # pGT875   RABALP1A
        alignment = alignments[5]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 42)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 42)
        self.assertEqual(alignment.shape, (2, 42 + 4))
        self.assertEqual(alignment.sequences[0].id, "RABALP1A")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 10)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.036)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 28.1)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [4973, 4987, 4987, 4990, 4990, 5002, 5003, 5016],
                        [240, 254, 256, 259, 260, 272, 272, 285],
                    ]
                ),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 5016)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "GCTGGAGAGAGCCA--TGG-TGGAGGCTGCGATGGAGGAGAGGATC")
        self.assertEqual(alignment[1], "GCCCGAAAGCACCACCTGGATGGAGAGACAGA-GGAGGAGAGGATC")
        # pGT875   OCDHPR
        alignment = alignments[6]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 12)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 12)
        self.assertEqual(alignment.shape, (2, 12 + 0))
        self.assertEqual(alignment.sequences[0].id, "OCDHPR")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 0)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.071)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 27.2)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[1499, 1511], [353, 341]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 1511)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "CCATGACCTGGT")
        self.assertEqual(alignment[1], "CCATGACCTGGT")
        # pGT875   RABALP1A
        alignment = alignments[7]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 12)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 12)
        self.assertEqual(alignment.shape, (2, 12 + 0))
        self.assertEqual(alignment.sequences[0].id, "RABALP1A")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 0)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.071)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 27.2)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[1499, 1511], [353, 341]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 1511)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "CCATGACCTGGT")
        self.assertEqual(alignment[1], "CCATGACCTGGT")
        # pGT875   RABGSTB
        alignment = alignments[8]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 16)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 16)
        self.assertEqual(alignment.shape, (2, 16 + 0))
        self.assertEqual(alignment.sequences[0].id, "RABGSTB")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 2)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.14)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 26.2)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[490, 506], [513, 497]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 506)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "CCTGGTTGAACTTCTC")
        self.assertEqual(alignment[1], "CCAGCTTGAACTTCTC")
        # pGT875   RABGLTR
        alignment = alignments[9]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 84)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 84)
        self.assertEqual(alignment.shape, (2, 84 + 0))
        self.assertEqual(alignment.sequences[0].id, "RABGLTR")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 39)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.28)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 25.2)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[1116, 1200], [499, 415]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 1200)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0],
            "GCATGGCTGGGTGGGGCAGGATTAGTGTGGGGGGAGTTGGGTGCTCAGGCAGGGCTATGAGGGATCTTGTTCATTTCCGGGCCC",
        )
        self.assertEqual(
            alignment[1],
            "GCAAGGTAGCGCAGGATGGCATTGCTCTGGGTGATCTTGTGTGATCCATCGATCAAGTAAGGCAGATTGGGAAAGTCCAGGCCC",
        )
        # pGT875   pGT875
        alignment = alignments[10]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 10)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 10)
        self.assertEqual(alignment.shape, (2, 10 + 0))
        self.assertEqual(alignment.sequences[0].id, "pGT875")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 0)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.28)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 25.2)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[792, 802], [357, 347]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 802)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "CCTGGTTCTC")
        self.assertEqual(alignment[1], "CCTGGTTCTC")
        # pGT875   BTGST
        alignment = alignments[11]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 71)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 71)
        self.assertEqual(alignment.shape, (2, 71 + 3))
        self.assertEqual(alignment.sequences[0].id, "BTGST")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 29)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.56)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 24.2)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [[280, 304, 304, 312, 312, 351], [353, 314, 312, 304, 303, 279]]
                ),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 351)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0],
            "CTGCGGCACCTGGGCCGCACCCTC-GGGCTGTA--TGGGAAGGACCAGCAGGAGGCGGCCCTGGTGGACATGGT",
        )
        self.assertEqual(
            alignment[1],
            "CTCTGGCTTCTGCTTCTCAAAGTCAGGGTTGTAACAGAGCATGATGAGCTGCATGCGGGTGTCCATGACCTGGT",
        )

    def test_m8CC(self):
        # Alignment file obtained by running
        # fasta36 -m 8CC seq/mgstm1.nt seq/gst.nlib
        # in the fasta36 source distribution
        path = "Fasta/nucleotide_m8CC.txt"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(
                alignments.commandline,
                "fasta36 -m 8CC seq/mgstm1.nt seq/gst.nlib",
            )
            alignments = list(alignments)
        self.assertEqual(len(alignments), 12)
        # pGT875   pGT875
        alignment = alignments[0]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 657)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 657)
        self.assertEqual(alignment.shape, (2, 657 + 0))
        self.assertEqual(alignment.sequences[0].id, "pGT875")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 0)
        self.assertAlmostEqual(alignment.annotations["evalue"], 3.6e-194)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 655.6)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[37, 694], [0, 657]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 694)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0],
            "ATGCCTATGATACTGGGATACTGGAACGTCCGCGGACTGACACACCCGATCCGCATGCTCCTGGAATACACAGACTCAAGCTATGATGAGAAGAGATACACCATGGGTGACGCTCCCGACTTTGACAGAAGCCAGTGGCTGAATGAGAAGTTCAAGCTGGGCCTGGACTTTCCCAATCTGCCTTACTTGATCGATGGATCACACAAGATCACCCAGAGCAATGCCATCCTGCGCTACCTTGCCCGAAAGCACCACCTGGATGGAGAGACAGAGGAGGAGAGGATCCGTGCAGACATTGTGGAGAACCAGGTCATGGACACCCGCATGCAGCTCATCATGCTCTGTTACAACCCTGACTTTGAGAAGCAGAAGCCAGAGTTCTTGAAGACCATCCCTGAGAAAATGAAGCTCTACTCTGAGTTCCTGGGCAAGAGGCCATGGTTTGCAGGGGACAAGGTCACCTATGTGGATTTCCTTGCTTATGACATTCTTGACCAGTACCGTATGTTTGAGCCCAAGTGCCTGGACGCCTTCCCAAACCTGAGGGACTTCCTGGCCCGCTTCGAGGGCCTCAAGAAGATCTCTGCCTACATGAAGAGTAGCCGCTACATCGCAACACCTATATTTTCAAAGATGGCCCACTGGAGTAACAAGTAG",
        )
        self.assertEqual(
            alignment[1],
            "ATGCCTATGATACTGGGATACTGGAACGTCCGCGGACTGACACACCCGATCCGCATGCTCCTGGAATACACAGACTCAAGCTATGATGAGAAGAGATACACCATGGGTGACGCTCCCGACTTTGACAGAAGCCAGTGGCTGAATGAGAAGTTCAAGCTGGGCCTGGACTTTCCCAATCTGCCTTACTTGATCGATGGATCACACAAGATCACCCAGAGCAATGCCATCCTGCGCTACCTTGCCCGAAAGCACCACCTGGATGGAGAGACAGAGGAGGAGAGGATCCGTGCAGACATTGTGGAGAACCAGGTCATGGACACCCGCATGCAGCTCATCATGCTCTGTTACAACCCTGACTTTGAGAAGCAGAAGCCAGAGTTCTTGAAGACCATCCCTGAGAAAATGAAGCTCTACTCTGAGTTCCTGGGCAAGAGGCCATGGTTTGCAGGGGACAAGGTCACCTATGTGGATTTCCTTGCTTATGACATTCTTGACCAGTACCGTATGTTTGAGCCCAAGTGCCTGGACGCCTTCCCAAACCTGAGGGACTTCCTGGCCCGCTTCGAGGGCCTCAAGAAGATCTCTGCCTACATGAAGAGTAGCCGCTACATCGCAACACCTATATTTTCAAAGATGGCCCACTGGAGTAACAAGTAG",
        )
        # pGT875   RABGLTR
        alignment = alignments[1]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 646)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 646)
        self.assertEqual(alignment.shape, (2, 646 + 0))
        self.assertEqual(alignment.sequences[0].id, "RABGLTR")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 135)
        self.assertAlmostEqual(alignment.annotations["evalue"], 1.9e-118)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 408.0)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[33, 679], [0, 646]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 679)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0],
            "ATGCCCATGACGCTGGGTTACTGGGACGTCCGTGGGCTGGCTCTGCCAATCCGCATGCTCCTGGAATACACGGACACCAGCTATGAGGAAAAGAAATACACCATGGGGGATGCTCCCAACTATGACCAAAGCAAGTGGCTGAGTGAGAAGTTCACCCTGGGCCTGGACTTTCCCAATCTGCCCTACCTAATTGATGGGACTCACAAGCTCACGCAGAGCAACGCCATCCTGCGCTACCTGGCCCGCAAGCACGGCCTGTGTGGGGAGACGGAAGAGGAGAGGATTCGCGTGGACATTCTGGAGAATCAGCTGATGGACAACCGCTTCCAACTTGTAAACGTCTGCTACAGTCCCGACTTTGAGAAGCTCAAGCCCGAGTACCTGAAGGGGCTCCCTGAGAAGCTGCAGCTGTACTCGCAGTTCCTGGGAAGCCTCCCCTGGTTCGCAGGGGACAAGATCACCTTCGCCGATTTCCTTGTCTACGACGTTCTTGACCAGAACCGGATATTTGTGCCTGGGTGCCTGGACGCGTTCCCAAACCTGAAGGACTTTCATGTCCGCTTTGAGGGCCTGCCGAAGATCTCTGCCTACATGAAGTCCAGCCGCTTTATCCGAGTCCCTGTGTTTTTAAAGAAGGCCACGTGGA",
        )
        self.assertEqual(
            alignment[1],
            "ATGCCTATGATACTGGGATACTGGAACGTCCGCGGACTGACACACCCGATCCGCATGCTCCTGGAATACACAGACTCAAGCTATGATGAGAAGAGATACACCATGGGTGACGCTCCCGACTTTGACAGAAGCCAGTGGCTGAATGAGAAGTTCAAGCTGGGCCTGGACTTTCCCAATCTGCCTTACTTGATCGATGGATCACACAAGATCACCCAGAGCAATGCCATCCTGCGCTACCTTGCCCGAAAGCACCACCTGGATGGAGAGACAGAGGAGGAGAGGATCCGTGCAGACATTGTGGAGAACCAGGTCATGGACACCCGCATGCAGCTCATCATGCTCTGTTACAACCCTGACTTTGAGAAGCAGAAGCCAGAGTTCTTGAAGACCATCCCTGAGAAAATGAAGCTCTACTCTGAGTTCCTGGGCAAGAGGCCATGGTTTGCAGGGGACAAGGTCACCTATGTGGATTTCCTTGCTTATGACATTCTTGACCAGTACCGTATGTTTGAGCCCAAGTGCCTGGACGCCTTCCCAAACCTGAGGGACTTCCTGGCCCGCTTCGAGGGCCTCAAGAAGATCTCTGCCTACATGAAGAGTAGCCGCTACATCGCAACACCTATATTTTCAAAGATGGCCCACTGGA",
        )
        # pGT875   BTGST
        alignment = alignments[2]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 413)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 413)
        self.assertEqual(alignment.shape, (2, 413 + 21))
        self.assertEqual(alignment.sequences[0].id, "BTGST")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 167)
        self.assertAlmostEqual(alignment.annotations["evalue"], 1.9e-07)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 45.7)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [
                            227,
                            376,
                            376,
                            383,
                            384,
                            401,
                            401,
                            461,
                            466,
                            472,
                            473,
                            486,
                            488,
                            501,
                            505,
                            535,
                            537,
                            543,
                            543,
                            655,
                        ],
                        [
                            175,
                            324,
                            325,
                            332,
                            332,
                            349,
                            352,
                            412,
                            412,
                            418,
                            418,
                            431,
                            431,
                            444,
                            444,
                            474,
                            474,
                            480,
                            482,
                            594,
                        ],
                    ]
                ),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 655)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0],
            "AGCTCCCCAAGTTCCAGGACGGAGACCTCACGCTGTACCAGTCCAATGCCATCCTGCGGCACCTGGGCCGCACCCTCGGGCTGTATGGGAAGGACCAGCAGGAGGCGGCCCTGGTGGACATGGTGAATGACGGTGTAGAGGACCTTCGC-TGCAAATACGTCTCCCTCATTTACA---CCAACTACGAGGCGGGCAAGGAGGACTATGTGAAGGCGCTGCCCCAGCACCTGAAGCCTTTCGAGACCCTGCTGTCCCAGAACAAGGGTGGCCAGGCCTTCATCGTGGGCGACCAGATCTCCTTTGCGGACTACAACCTGCT--GGACCTGCTTCGGATTCACCAGGTCCTGGCCCCCAGCTGTCTGGACTCCTTCCCCCTGCTCTCAGCCTACGTGGCCCGTCTCAACTCCCGGCCCAAGCTCAAGGCCTTCCTG",
        )
        self.assertEqual(
            alignment[1],
            "ATCTGCCTTACTTGATCGATGGATCACACAAGATCACCCAGAGCAATGCCATCCTGCGCTACCTTGCCCGAAAGCACCACCTGGATGGAGAGACAGAGGAGGAGAGGATCCGTGCAGACATTGTGGAGAACCAGGTCATGGACACCCGCATGCAGCT-CATCATGCTCTGTTACAACCCTGACTTTGAGAAGCAGAAGCCAGAGTTCTTGAAGACCATCCCTGAGAAAATGAAGCTCT-----ACTCTG-AGTTCCTGGGCAA--GAGGCCATGGTTT----GCAGGGGACAAGGTCACCTATGTGGATTTC--CTTGCTTATGACATTCTTGACCAGTACCGTATGTTTGAGCCCAAGTGCCTGGACGCCTTCCCAAACCTGAGGGACTTCCTGGCCCGCTTCGAGGGCCTCAAGAAGATCTCTGCCTACATG",
        )
        # pGT875   RABGSTB
        alignment = alignments[3]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 127)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 127)
        self.assertEqual(alignment.shape, (2, 127 + 8))
        self.assertEqual(alignment.sequences[0].id, "RABGSTB")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 42)
        self.assertAlmostEqual(alignment.annotations["evalue"], 3.2e-07)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 45.0)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [
                            156,
                            171,
                            173,
                            190,
                            190,
                            201,
                            202,
                            260,
                            261,
                            272,
                            272,
                            279,
                            279,
                            287,
                        ],
                        [
                            158,
                            173,
                            173,
                            190,
                            192,
                            203,
                            203,
                            261,
                            261,
                            272,
                            273,
                            280,
                            281,
                            289,
                        ],
                    ]
                ),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 287)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0],
            "GGGTATTGATGTTCCAGCAAGTGCCCATGGTTGA--GATTGATGGGATGAAGCTGGTGCAGACCAGAGCCATTTTCAACTACATTGCAGACAAGCACAACCTGTATGGGAAAGACATA-AAGGAGA-GAGCCCTG",
        )
        self.assertEqual(
            alignment[1],
            "GGGCCTGGACTTTCC--CAATCTGCCTTACTTGATCGATGGATCACA-CAAGATCACCCAGAGCAATGCCATCCTGCGCTACCTTGCCCGAAAGCACCACCTGGAT-GGAGAGACAGAGGAGGAGAGGATCCGTG",
        )
        # pGT875   OCDHPR
        alignment = alignments[4]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 23)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 23)
        self.assertEqual(alignment.shape, (2, 23 + 1))
        self.assertEqual(alignment.sequences[0].id, "OCDHPR")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 2)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.012)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 29.7)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[2302, 2319, 2319, 2325], [265, 282, 283, 289]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 2325)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "AGACAGAGGAGGAGAAG-TCTGTG")
        self.assertEqual(alignment[1], "AGACAGAGGAGGAGAGGATCCGTG")
        # pGT875   RABALP1A
        alignment = alignments[5]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 42)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 42)
        self.assertEqual(alignment.shape, (2, 42 + 4))
        self.assertEqual(alignment.sequences[0].id, "RABALP1A")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 10)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.046)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 27.8)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [
                        [4973, 4987, 4987, 4990, 4990, 5002, 5003, 5016],
                        [240, 254, 256, 259, 260, 272, 272, 285],
                    ]
                ),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 5016)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "GCTGGAGAGAGCCA--TGG-TGGAGGCTGCGATGGAGGAGAGGATC")
        self.assertEqual(alignment[1], "GCCCGAAAGCACCACCTGGATGGAGAGACAGA-GGAGGAGAGGATC")
        # pGT875   OCDHPR
        alignment = alignments[6]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 12)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 12)
        self.assertEqual(alignment.shape, (2, 12 + 0))
        self.assertEqual(alignment.sequences[0].id, "OCDHPR")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 0)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.09)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 26.8)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[1499, 1511], [353, 341]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 1511)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "CCATGACCTGGT")
        self.assertEqual(alignment[1], "CCATGACCTGGT")
        # pGT875   RABALP1A
        alignment = alignments[7]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 12)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 12)
        self.assertEqual(alignment.shape, (2, 12 + 0))
        self.assertEqual(alignment.sequences[0].id, "RABALP1A")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 0)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.09)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 26.8)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[1499, 1511], [353, 341]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 1511)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "CCATGACCTGGT")
        self.assertEqual(alignment[1], "CCATGACCTGGT")
        # pGT875   RABGSTB
        alignment = alignments[8]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 16)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 16)
        self.assertEqual(alignment.shape, (2, 16 + 0))
        self.assertEqual(alignment.sequences[0].id, "RABGSTB")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 2)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.18)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 25.8)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[490, 506], [513, 497]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 506)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "CCTGGTTGAACTTCTC")
        self.assertEqual(alignment[1], "CCAGCTTGAACTTCTC")
        # pGT875   RABGLTR
        alignment = alignments[9]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 84)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 84)
        self.assertEqual(alignment.shape, (2, 84 + 0))
        self.assertEqual(alignment.sequences[0].id, "RABGLTR")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 39)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.35)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 24.9)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[1116, 1200], [499, 415]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 1200)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0],
            "GCATGGCTGGGTGGGGCAGGATTAGTGTGGGGGGAGTTGGGTGCTCAGGCAGGGCTATGAGGGATCTTGTTCATTTCCGGGCCC",
        )
        self.assertEqual(
            alignment[1],
            "GCAAGGTAGCGCAGGATGGCATTGCTCTGGGTGATCTTGTGTGATCCATCGATCAAGTAAGGCAGATTGGGAAAGTCCAGGCCC",
        )
        # pGT875   pGT875
        alignment = alignments[10]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 10)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 10)
        self.assertEqual(alignment.shape, (2, 10 + 0))
        self.assertEqual(alignment.sequences[0].id, "pGT875")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 0)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.35)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 24.9)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[792, 802], [357, 347]]),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 802)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(alignment[0], "CCTGGTTCTC")
        self.assertEqual(alignment[1], "CCTGGTTCTC")
        # pGT875   BTGST
        alignment = alignments[11]
        self.assertEqual(alignment.annotations["program"], "FASTA 36.3.8h May, 2020")
        self.assertEqual(alignment.annotations["database"], "seq/gst.nlib")
        self.assertEqual(len(alignment), 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 71)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 71)
        self.assertEqual(alignment.shape, (2, 71 + 3))
        self.assertEqual(alignment.sequences[0].id, "BTGST")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 29)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.68)
        self.assertAlmostEqual(alignment.annotations["bit_score"], 23.9)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [[280, 304, 304, 312, 312, 351], [353, 314, 312, 304, 303, 279]]
                ),
            )
        )
        query = self.query
        target = self.targets[alignment.sequences[0].id]
        self.assertEqual(len(alignment.sequences[0].seq), 351)
        self.assertEqual(len(alignment.sequences[1].seq), len(query))
        alignment.sequences[0].seq = target
        alignment.sequences[1].seq = query
        self.assertEqual(
            alignment[0],
            "CTGCGGCACCTGGGCCGCACCCTC-GGGCTGTA--TGGGAAGGACCAGCAGGAGGCGGCCCTGGTGGACATGGT",
        )
        self.assertEqual(
            alignment[1],
            "CTCTGGCTTCTGCTTCTCAAAGTCAGGGTTGTAACAGAGCATGATGAGCTGCATGCGGGTGTCCATGACCTGGT",
        )


class TestFastaBasic(unittest.TestCase):
    def test_empty(self):
        import io

        stream = io.StringIO()
        with self.assertRaisesRegex(ValueError, "Empty file."):
            AlignmentIterator(stream)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
