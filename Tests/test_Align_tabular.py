# Copyright 2021 by Michiel de Hoon.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Bio.Align.tabular module."""
import unittest
import warnings
import os

from Bio.Seq import Seq
from Bio import SeqIO
from Bio import BiopythonExperimentalWarning

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    from Bio.Align.tabular import AlignmentIterator


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
            self.assertEqual(alignments.program, "FASTA")
            self.assertEqual(alignments.version, "36.3.8h May, 2020")
            self.assertEqual(alignments.database, "seq/prot_test.lseg")
            alignments = list(alignments)
        self.assertEqual(len(alignments), 12)
        # sp|P10649|GSTM1_MOUSE   sp|P09488|GSTM1_HUMAN
        alignment = alignments[0]
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 77.98)
        self.assertEqual(alignment.annotations["gap opens"], 0)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 218)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 218)
        self.assertEqual(alignment.shape, (2, 218 + 0))
        self.assertEqual(alignment.sequences[0].id, "sp|P09488|GSTM1_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 48)
        self.assertAlmostEqual(alignment.annotations["evalue"], 6.1e-78)
        self.assertAlmostEqual(alignment.annotations["bit score"], 275.6)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 29.76)
        self.assertEqual(alignment.annotations["gap opens"], 18)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 205)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 205)
        self.assertEqual(alignment.shape, (2, 205 + 18))
        self.assertEqual(alignment.sequences[0].id, "sp|P00502|GSTA1_RAT")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 144)
        self.assertAlmostEqual(alignment.annotations["evalue"], 3.1e-13)
        self.assertAlmostEqual(alignment.annotations["bit score"], 60.7)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 27.03)
        self.assertEqual(alignment.annotations["gap opens"], 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 37)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 37)
        self.assertEqual(alignment.shape, (2, 37 + 2))
        self.assertEqual(alignment.sequences[0].id, "sp|P69905|HBA_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 27)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.19)
        self.assertAlmostEqual(alignment.annotations["bit score"], 20.9)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 24.29)
        self.assertEqual(alignment.annotations["gap opens"], 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 70)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 70)
        self.assertEqual(alignment.shape, (2, 70 + 2))
        self.assertEqual(alignment.sequences[0].id, "sp|P00517|KAPCA_BOVIN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 53)
        self.assertAlmostEqual(alignment.annotations["evalue"], 1.4)
        self.assertAlmostEqual(alignment.annotations["bit score"], 19.2)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 57.14)
        self.assertEqual(alignment.annotations["gap opens"], 0)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 7)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 7)
        self.assertEqual(alignment.shape, (2, 7 + 0))
        self.assertEqual(alignment.sequences[0].id, "sp|P14960|RBS_GUITH")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 3)
        self.assertAlmostEqual(alignment.annotations["evalue"], 1.6)
        self.assertAlmostEqual(alignment.annotations["bit score"], 17.7)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 45.00)
        self.assertEqual(alignment.annotations["gap opens"], 10)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 40)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 40)
        self.assertEqual(alignment.shape, (2, 40 + 10))
        self.assertEqual(alignment.sequences[0].id, "sp|P01593|KV101_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 22)
        self.assertAlmostEqual(alignment.annotations["evalue"], 2.5)
        self.assertAlmostEqual(alignment.annotations["bit score"], 16.6)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 31.03)
        self.assertEqual(alignment.annotations["gap opens"], 10)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 58)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 58)
        self.assertEqual(alignment.shape, (2, 58 + 10))
        self.assertEqual(alignment.sequences[0].id, "sp|P99998|CYC_PANTR")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 40)
        self.assertAlmostEqual(alignment.annotations["evalue"], 2.7)
        self.assertAlmostEqual(alignment.annotations["bit score"], 16.4)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 31.11)
        self.assertEqual(alignment.annotations["gap opens"], 9)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 45)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 45)
        self.assertEqual(alignment.shape, (2, 45 + 9))
        self.assertEqual(alignment.sequences[0].id, "sp|P02585|TNNC2_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 31)
        self.assertAlmostEqual(alignment.annotations["evalue"], 3.9)
        self.assertAlmostEqual(alignment.annotations["bit score"], 16.4)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 66.67)
        self.assertEqual(alignment.annotations["gap opens"], 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 9)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 9)
        self.assertEqual(alignment.shape, (2, 9 + 2))
        self.assertEqual(alignment.sequences[0].id, "sp|P60615|NXL1A_BUNMU")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 3)
        self.assertAlmostEqual(alignment.annotations["evalue"], 4.2)
        self.assertAlmostEqual(alignment.annotations["bit score"], 15.6)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 50.00)
        self.assertEqual(alignment.annotations["gap opens"], 0)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 4)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 4)
        self.assertEqual(alignment.shape, (2, 4 + 0))
        self.assertEqual(alignment.sequences[0].id, "sp|P00193|FER_PEPAS")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 2)
        self.assertAlmostEqual(alignment.annotations["evalue"], 4.2)
        self.assertAlmostEqual(alignment.annotations["bit score"], 14.7)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 28.95)
        self.assertEqual(alignment.annotations["gap opens"], 1)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 38)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 38)
        self.assertEqual(alignment.shape, (2, 38 + 1))
        self.assertEqual(alignment.sequences[0].id, "sp|P03435|HEMA_I75A3")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 27)
        self.assertAlmostEqual(alignment.annotations["evalue"], 4.7)
        self.assertAlmostEqual(alignment.annotations["bit score"], 17.9)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 19.40)
        self.assertEqual(alignment.annotations["gap opens"], 4)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 67)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 67)
        self.assertEqual(alignment.shape, (2, 67 + 4))
        self.assertEqual(alignment.sequences[0].id, "sp|P01834|IGKC_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 54)
        self.assertAlmostEqual(alignment.annotations["evalue"], 6.3)
        self.assertAlmostEqual(alignment.annotations["bit score"], 14.9)
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
            self.assertEqual(alignments.program, "FASTA")
            self.assertEqual(alignments.version, "36.3.8h May, 2020")
            self.assertEqual(alignments.database, "seq/prot_test.lseg")
            alignments = list(alignments)
        self.assertEqual(len(alignments), 12)
        # sp|P10649|GSTM1_MOUSE   sp|P09488|GSTM1_HUMAN
        alignment = alignments[0]
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 77.98)
        self.assertEqual(alignment.annotations["gap opens"], 0)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 218)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 218)
        self.assertEqual(alignment.shape, (2, 218 + 0))
        self.assertEqual(alignment.sequences[0].id, "sp|P09488|GSTM1_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 48)
        self.assertAlmostEqual(alignment.annotations["evalue"], 7.6e-83)
        self.assertAlmostEqual(alignment.annotations["bit score"], 291.9)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 29.76)
        self.assertEqual(alignment.annotations["gap opens"], 18)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 205)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 205)
        self.assertEqual(alignment.shape, (2, 205 + 18))
        self.assertEqual(alignment.sequences[0].id, "sp|P00502|GSTA1_RAT")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 144)
        self.assertAlmostEqual(alignment.annotations["evalue"], 4.4e-14)
        self.assertAlmostEqual(alignment.annotations["bit score"], 63.5)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 27.03)
        self.assertEqual(alignment.annotations["gap opens"], 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 37)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 37)
        self.assertEqual(alignment.shape, (2, 37 + 2))
        self.assertEqual(alignment.sequences[0].id, "sp|P69905|HBA_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 27)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.15)
        self.assertAlmostEqual(alignment.annotations["bit score"], 21.2)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 24.29)
        self.assertEqual(alignment.annotations["gap opens"], 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 70)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 70)
        self.assertEqual(alignment.shape, (2, 70 + 2))
        self.assertEqual(alignment.sequences[0].id, "sp|P00517|KAPCA_BOVIN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 53)
        self.assertAlmostEqual(alignment.annotations["evalue"], 1.2)
        self.assertAlmostEqual(alignment.annotations["bit score"], 19.4)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 57.14)
        self.assertEqual(alignment.annotations["gap opens"], 0)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 7)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 7)
        self.assertEqual(alignment.shape, (2, 7 + 0))
        self.assertEqual(alignment.sequences[0].id, "sp|P14960|RBS_GUITH")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 3)
        self.assertAlmostEqual(alignment.annotations["evalue"], 1.5)
        self.assertAlmostEqual(alignment.annotations["bit score"], 17.8)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 45.00)
        self.assertEqual(alignment.annotations["gap opens"], 10)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 40)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 40)
        self.assertEqual(alignment.shape, (2, 40 + 10))
        self.assertEqual(alignment.sequences[0].id, "sp|P01593|KV101_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 22)
        self.assertAlmostEqual(alignment.annotations["evalue"], 2.4)
        self.assertAlmostEqual(alignment.annotations["bit score"], 16.7)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 31.03)
        self.assertEqual(alignment.annotations["gap opens"], 10)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 58)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 58)
        self.assertEqual(alignment.shape, (2, 58 + 10))
        self.assertEqual(alignment.sequences[0].id, "sp|P99998|CYC_PANTR")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 40)
        self.assertAlmostEqual(alignment.annotations["evalue"], 2.7)
        self.assertAlmostEqual(alignment.annotations["bit score"], 16.5)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 31.11)
        self.assertEqual(alignment.annotations["gap opens"], 9)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 45)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 45)
        self.assertEqual(alignment.shape, (2, 45 + 9))
        self.assertEqual(alignment.sequences[0].id, "sp|P02585|TNNC2_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 31)
        self.assertAlmostEqual(alignment.annotations["evalue"], 3.8)
        self.assertAlmostEqual(alignment.annotations["bit score"], 16.5)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 66.67)
        self.assertEqual(alignment.annotations["gap opens"], 2)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 9)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 9)
        self.assertEqual(alignment.shape, (2, 9 + 2))
        self.assertEqual(alignment.sequences[0].id, "sp|P60615|NXL1A_BUNMU")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 3)
        self.assertAlmostEqual(alignment.annotations["evalue"], 4.2)
        self.assertAlmostEqual(alignment.annotations["bit score"], 15.6)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 28.95)
        self.assertEqual(alignment.annotations["gap opens"], 1)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 38)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 38)
        self.assertEqual(alignment.shape, (2, 38 + 1))
        self.assertEqual(alignment.sequences[0].id, "sp|P03435|HEMA_I75A3")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 27)
        self.assertAlmostEqual(alignment.annotations["evalue"], 4.4)
        self.assertAlmostEqual(alignment.annotations["bit score"], 18.1)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 50.00)
        self.assertEqual(alignment.annotations["gap opens"], 0)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 4)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 4)
        self.assertEqual(alignment.shape, (2, 4 + 0))
        self.assertEqual(alignment.sequences[0].id, "sp|P00193|FER_PEPAS")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 2)
        self.assertAlmostEqual(alignment.annotations["evalue"], 4.4)
        self.assertAlmostEqual(alignment.annotations["bit score"], 14.7)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 19.40)
        self.assertEqual(alignment.annotations["gap opens"], 4)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 67)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 67)
        self.assertEqual(alignment.shape, (2, 67 + 4))
        self.assertEqual(alignment.sequences[0].id, "sp|P01834|IGKC_HUMAN")
        self.assertEqual(alignment.sequences[1].id, "sp|P10649|GSTM1_MOUSE")
        self.assertEqual(alignment.annotations["mismatches"], 54)
        self.assertAlmostEqual(alignment.annotations["evalue"], 6.4)
        self.assertAlmostEqual(alignment.annotations["bit score"], 14.9)
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
            self.assertEqual(alignments.program, "FASTA")
            self.assertEqual(alignments.version, "36.3.8h May, 2020")
            self.assertEqual(alignments.database, "seq/gst.nlib")
            alignments = list(alignments)
        self.assertEqual(len(alignments), 12)
        # pGT875   pGT875
        alignment = alignments[0]
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 100.00)
        self.assertEqual(alignment.annotations["gap opens"], 0)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 657)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 657)
        self.assertEqual(alignment.shape, (2, 657 + 0))
        self.assertEqual(alignment.sequences[0].id, "pGT875")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 0)
        self.assertAlmostEqual(alignment.annotations["evalue"], 3.6e-194)
        self.assertAlmostEqual(alignment.annotations["bit score"], 666.0)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 79.10)
        self.assertEqual(alignment.annotations["gap opens"], 0)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 646)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 646)
        self.assertEqual(alignment.shape, (2, 646 + 0))
        self.assertEqual(alignment.sequences[0].id, "RABGLTR")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 135)
        self.assertAlmostEqual(alignment.annotations["evalue"], 1.9e-118)
        self.assertAlmostEqual(alignment.annotations["bit score"], 414.4)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 59.56)
        self.assertEqual(alignment.annotations["gap opens"], 21)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 413)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 413)
        self.assertEqual(alignment.shape, (2, 413 + 21))
        self.assertEqual(alignment.sequences[0].id, "BTGST")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 167)
        self.assertAlmostEqual(alignment.annotations["evalue"], 1.2e-07)
        self.assertAlmostEqual(alignment.annotations["bit score"], 46.4)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 66.93)
        self.assertEqual(alignment.annotations["gap opens"], 8)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 127)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 127)
        self.assertEqual(alignment.shape, (2, 127 + 8))
        self.assertEqual(alignment.sequences[0].id, "RABGSTB")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 42)
        self.assertAlmostEqual(alignment.annotations["evalue"], 2.1e-07)
        self.assertAlmostEqual(alignment.annotations["bit score"], 45.6)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 91.30)
        self.assertEqual(alignment.annotations["gap opens"], 1)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 23)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 23)
        self.assertEqual(alignment.shape, (2, 23 + 1))
        self.assertEqual(alignment.sequences[0].id, "OCDHPR")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 2)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.0092)
        self.assertAlmostEqual(alignment.annotations["bit score"], 30.1)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 76.19)
        self.assertEqual(alignment.annotations["gap opens"], 4)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 42)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 42)
        self.assertEqual(alignment.shape, (2, 42 + 4))
        self.assertEqual(alignment.sequences[0].id, "RABALP1A")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 10)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.036)
        self.assertAlmostEqual(alignment.annotations["bit score"], 28.1)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 100.00)
        self.assertEqual(alignment.annotations["gap opens"], 0)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 12)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 12)
        self.assertEqual(alignment.shape, (2, 12 + 0))
        self.assertEqual(alignment.sequences[0].id, "OCDHPR")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 0)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.071)
        self.assertAlmostEqual(alignment.annotations["bit score"], 27.2)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[1499, 1511], [316, 304]]),
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 100.00)
        self.assertEqual(alignment.annotations["gap opens"], 0)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 12)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 12)
        self.assertEqual(alignment.shape, (2, 12 + 0))
        self.assertEqual(alignment.sequences[0].id, "RABALP1A")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 0)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.071)
        self.assertAlmostEqual(alignment.annotations["bit score"], 27.2)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[1499, 1511], [316, 304]]),
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 87.50)
        self.assertEqual(alignment.annotations["gap opens"], 0)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 16)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 16)
        self.assertEqual(alignment.shape, (2, 16 + 0))
        self.assertEqual(alignment.sequences[0].id, "RABGSTB")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 2)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.14)
        self.assertAlmostEqual(alignment.annotations["bit score"], 26.2)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[490, 506], [160, 144]]),
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 53.57)
        self.assertEqual(alignment.annotations["gap opens"], 0)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 84)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 84)
        self.assertEqual(alignment.shape, (2, 84 + 0))
        self.assertEqual(alignment.sequences[0].id, "RABGLTR")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 39)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.28)
        self.assertAlmostEqual(alignment.annotations["bit score"], 25.2)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[1116, 1200], [242, 158]]),
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 100.00)
        self.assertEqual(alignment.annotations["gap opens"], 0)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 10)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 10)
        self.assertEqual(alignment.shape, (2, 10 + 0))
        self.assertEqual(alignment.sequences[0].id, "pGT875")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 0)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.28)
        self.assertAlmostEqual(alignment.annotations["bit score"], 25.2)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[792, 802], [310, 300]]),
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 59.15)
        self.assertEqual(alignment.annotations["gap opens"], 3)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 71)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 71)
        self.assertEqual(alignment.shape, (2, 71 + 3))
        self.assertEqual(alignment.sequences[0].id, "BTGST")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 29)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.56)
        self.assertAlmostEqual(alignment.annotations["bit score"], 24.2)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [[280, 304, 304, 312, 312, 351], [378, 354, 353, 345, 343, 304]]
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
            self.assertEqual(alignments.program, "FASTA")
            self.assertEqual(alignments.version, "36.3.8h May, 2020")
            self.assertEqual(alignments.database, "seq/gst.nlib")
            alignments = list(alignments)
        self.assertEqual(len(alignments), 12)
        # pGT875   pGT875
        alignment = alignments[0]
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 100.00)
        self.assertEqual(alignment.annotations["gap opens"], 0)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 657)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 657)
        self.assertEqual(alignment.shape, (2, 657 + 0))
        self.assertEqual(alignment.sequences[0].id, "pGT875")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 0)
        self.assertAlmostEqual(alignment.annotations["evalue"], 3.6e-194)
        self.assertAlmostEqual(alignment.annotations["bit score"], 655.6)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 79.10)
        self.assertEqual(alignment.annotations["gap opens"], 0)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 646)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 646)
        self.assertEqual(alignment.shape, (2, 646 + 0))
        self.assertEqual(alignment.sequences[0].id, "RABGLTR")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 135)
        self.assertAlmostEqual(alignment.annotations["evalue"], 1.9e-118)
        self.assertAlmostEqual(alignment.annotations["bit score"], 408.0)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 59.56)
        self.assertEqual(alignment.annotations["gap opens"], 21)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 413)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 413)
        self.assertEqual(alignment.shape, (2, 413 + 21))
        self.assertEqual(alignment.sequences[0].id, "BTGST")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 167)
        self.assertAlmostEqual(alignment.annotations["evalue"], 1.9e-07)
        self.assertAlmostEqual(alignment.annotations["bit score"], 45.7)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 66.93)
        self.assertEqual(alignment.annotations["gap opens"], 8)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 127)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 127)
        self.assertEqual(alignment.shape, (2, 127 + 8))
        self.assertEqual(alignment.sequences[0].id, "RABGSTB")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 42)
        self.assertAlmostEqual(alignment.annotations["evalue"], 3.2e-07)
        self.assertAlmostEqual(alignment.annotations["bit score"], 45.0)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 91.30)
        self.assertEqual(alignment.annotations["gap opens"], 1)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 23)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 23)
        self.assertEqual(alignment.shape, (2, 23 + 1))
        self.assertEqual(alignment.sequences[0].id, "OCDHPR")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 2)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.012)
        self.assertAlmostEqual(alignment.annotations["bit score"], 29.7)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 76.19)
        self.assertEqual(alignment.annotations["gap opens"], 4)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 42)
        self.assertEqual(sum(end - start for start, end in aligned[1]), 42)
        self.assertEqual(alignment.shape, (2, 42 + 4))
        self.assertEqual(alignment.sequences[0].id, "RABALP1A")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 10)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.046)
        self.assertAlmostEqual(alignment.annotations["bit score"], 27.8)
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 100.00)
        self.assertEqual(alignment.annotations["gap opens"], 0)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 12)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 12)
        self.assertEqual(alignment.shape, (2, 12 + 0))
        self.assertEqual(alignment.sequences[0].id, "OCDHPR")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 0)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.09)
        self.assertAlmostEqual(alignment.annotations["bit score"], 26.8)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[1499, 1511], [316, 304]]),
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 100.00)
        self.assertEqual(alignment.annotations["gap opens"], 0)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 12)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 12)
        self.assertEqual(alignment.shape, (2, 12 + 0))
        self.assertEqual(alignment.sequences[0].id, "RABALP1A")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 0)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.09)
        self.assertAlmostEqual(alignment.annotations["bit score"], 26.8)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[1499, 1511], [316, 304]]),
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 87.50)
        self.assertEqual(alignment.annotations["gap opens"], 0)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 16)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 16)
        self.assertEqual(alignment.shape, (2, 16 + 0))
        self.assertEqual(alignment.sequences[0].id, "RABGSTB")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 2)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.18)
        self.assertAlmostEqual(alignment.annotations["bit score"], 25.8)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[490, 506], [160, 144]]),
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 53.57)
        self.assertEqual(alignment.annotations["gap opens"], 0)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 84)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 84)
        self.assertEqual(alignment.shape, (2, 84 + 0))
        self.assertEqual(alignment.sequences[0].id, "RABGLTR")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 39)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.35)
        self.assertAlmostEqual(alignment.annotations["bit score"], 24.9)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[1116, 1200], [242, 158]]),
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 100.00)
        self.assertEqual(alignment.annotations["gap opens"], 0)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 10)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 10)
        self.assertEqual(alignment.shape, (2, 10 + 0))
        self.assertEqual(alignment.sequences[0].id, "pGT875")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 0)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.35)
        self.assertAlmostEqual(alignment.annotations["bit score"], 24.9)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[792, 802], [310, 300]]),
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
        self.assertEqual(len(alignment), 2)
        self.assertAlmostEqual(alignment.annotations["% identity"], 59.15)
        self.assertEqual(alignment.annotations["gap opens"], 3)
        aligned = alignment.aligned
        self.assertEqual(sum(end - start for start, end in aligned[0]), 71)
        self.assertEqual(sum(start - end for start, end in aligned[1]), 71)
        self.assertEqual(alignment.shape, (2, 71 + 3))
        self.assertEqual(alignment.sequences[0].id, "BTGST")
        self.assertEqual(alignment.sequences[1].id, "pGT875")
        self.assertEqual(alignment.annotations["mismatches"], 29)
        self.assertAlmostEqual(alignment.annotations["evalue"], 0.68)
        self.assertAlmostEqual(alignment.annotations["bit score"], 23.9)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array(
                    [[280, 304, 304, 312, 312, 351], [378, 354, 353, 345, 343, 304]]
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


class TestBlast(unittest.TestCase):

    def test_tblastn_011(self):
        path = "Blast/tab_2226_tblastn_011.txt"
        with open(path) as stream:
            alignments = AlignmentIterator(stream)
            self.assertEqual(alignments.program, "TBLASTN")
            self.assertEqual(alignments.version, "2.2.26+")
            alignment = next(alignments)
            self.assertEqual(alignment.query.id, "gi|16080617|ref|NP_391444.1|")
            self.assertEqual(alignment.query.annotations["gi"], "0")
            self.assertEqual(alignment.query.annotations["acc."], "gi|16080617|ref|NP_391444.1|")
            self.assertEqual(alignment.query.annotations["acc.ver"], "gi|16080617|ref|NP_391444.1|")
            self.assertEqual(len(alignment.query.seq), 102)
            self.assertEqual(alignment.target.id, "gi|145479850|ref|XM_001425911.1|")
            self.assertEqual(alignment.target.annotations["ids"], "gi|145479850|ref|XM_001425911.1|")
            self.assertEqual(alignment.target.annotations["gi"], "0")
            self.assertEqual(alignment.target.annotations["gis"], "0")
            self.assertEqual(alignment.target.annotations["acc."], "gi|145479850|ref|XM_001425911.1|")
            self.assertEqual(alignment.target.annotations["acc.ver"], "gi|145479850|ref|XM_001425911.1|")
            self.assertEqual(alignment.target.annotations["length"], 4632)
            self.assertTrue(
                numpy.array_equal(
                    alignment.coordinates, numpy.array([[0, 43], [30, 73]]),
                )
            )
            self.assertEqual(alignment.target.annotations["start"], 1743)
            self.assertEqual(alignment.target.annotations["end"], 1872)
            self.assertEqual(alignment.target.seq, "PKTATGTKKGTIIGLLSIHTILFILTSHALSLEVKEQT*KDID")
            self.assertEqual(alignment.query.seq[30:73], "PDSNIETKEGTYVGLADTHTIEVTVDNEPVSLDITEESTSDLD")
            self.assertAlmostEqual(alignment.annotations["evalue"], 1e-05)
            self.assertAlmostEqual(alignment.annotations["bit score"], 34.7)
            self.assertAlmostEqual(alignment.score, 78)
            self.assertEqual(alignment.shape, (2, 43))
            self.assertAlmostEqual(alignment.annotations["% identity"], 34.88)
            self.assertEqual(alignment.annotations["identical"], 15)
            self.assertEqual(alignment.annotations["mismatches"], 28)
            self.assertEqual(alignment.annotations["positives"], 26)
            self.assertEqual(alignment.annotations["gap opens"], 0)
            self.assertEqual(alignment.annotations["gaps"], 0)
            self.assertAlmostEqual(alignment.annotations["% positives"], 60.47)
            self.assertEqual(alignment.annotations["query/sbjct frames"], "0/1")
            self.assertEqual(alignment.query.annotations["frame"], "0")
            self.assertEqual(alignment.target.annotations["frame"], "1")
            alignment = next(alignments)
            self.assertEqual(alignment.query.id, "gi|16080617|ref|NP_391444.1|")
            self.assertEqual(alignment.query.annotations["gi"], "0")
            self.assertEqual(alignment.query.annotations["acc."], "gi|16080617|ref|NP_391444.1|")
            self.assertEqual(alignment.query.annotations["acc.ver"], "gi|16080617|ref|NP_391444.1|")
            self.assertEqual(len(alignment.query.seq), 102)
            self.assertEqual(alignment.target.id, "gi|72012412|ref|XM_777959.1|")
            self.assertEqual(alignment.target.annotations["ids"], "gi|72012412|ref|XM_777959.1|")
            self.assertEqual(alignment.target.annotations["gi"], "0")
            self.assertEqual(alignment.target.annotations["gis"], "0")
            self.assertEqual(alignment.target.annotations["acc."], "gi|72012412|ref|XM_777959.1|")
            self.assertEqual(alignment.target.annotations["acc.ver"], "gi|72012412|ref|XM_777959.1|")
            self.assertEqual(alignment.target.annotations["length"], 1593)
            self.assertTrue(
                numpy.array_equal(
                    alignment.coordinates,
                    numpy.array([[ 0, 35, 43, 59], [43, 78, 78, 94]])
                )
            )
            self.assertEqual(alignment.target.annotations["start"], 1056)
            self.assertEqual(alignment.target.annotations["end"], 1233)
            self.assertEqual(alignment.target.seq, "GLVPDHTLILPVGHYQSMLDLTEEVQTELDQFKSALRKYYLSKGKTCVIYERNFRTQHL")
            self.assertEqual(alignment.query.seq[43:94], "GLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLL")
            self.assertAlmostEqual(alignment.annotations["evalue"], 1e-04)
            self.assertAlmostEqual(alignment.annotations["bit score"], 31.6)
            self.assertAlmostEqual(alignment.score, 70)
            self.assertEqual(alignment.shape, (2, 59))
            self.assertAlmostEqual(alignment.annotations["% identity"], 33.90)
            self.assertEqual(alignment.annotations["identical"], 20)
            self.assertEqual(alignment.annotations["mismatches"], 31)
            self.assertEqual(alignment.annotations["positives"], 29)
            self.assertEqual(alignment.annotations["gap opens"], 1)
            self.assertEqual(alignment.annotations["gaps"], 8)
            self.assertAlmostEqual(alignment.annotations["% positives"], 49.15)
            self.assertEqual(alignment.annotations["query/sbjct frames"], "0/1")
            self.assertEqual(alignment.query.annotations["frame"], "0")
            self.assertEqual(alignment.target.annotations["frame"], "1")
            alignment = next(alignments)
            self.assertEqual(alignment.query.id, "gi|16080617|ref|NP_391444.1|")
            self.assertEqual(alignment.query.annotations["gi"], "0")
            self.assertEqual(alignment.query.annotations["acc."], "gi|16080617|ref|NP_391444.1|")
            self.assertEqual(alignment.query.annotations["acc.ver"], "gi|16080617|ref|NP_391444.1|")
            self.assertEqual(len(alignment.query.seq), 102)
            self.assertEqual(alignment.target.id, "gi|115975252|ref|XM_001180111.1|")
            self.assertEqual(alignment.target.annotations["ids"], "gi|115975252|ref|XM_001180111.1|")
            self.assertEqual(alignment.target.annotations["gi"], "0")
            self.assertEqual(alignment.target.annotations["gis"], "0")
            self.assertEqual(alignment.target.annotations["acc."], "gi|115975252|ref|XM_001180111.1|")
            self.assertEqual(alignment.target.annotations["acc.ver"], "gi|115975252|ref|XM_001180111.1|")
            self.assertEqual(alignment.target.annotations["length"], 1593)
            self.assertTrue(
                numpy.array_equal(
                    alignment.coordinates,
                    numpy.array([[0, 35, 43, 59], [43, 78, 78, 94]])
                )
            )
            self.assertEqual(alignment.target.annotations["start"], 1056)
            self.assertEqual(alignment.target.annotations["end"], 1233)
            self.assertEqual(alignment.target.seq, "GLVPDHTLILPVGHYQSMLDLTEEVQTELDQFKSALRKYYLSKGKTCVIYERNFRTQHL")
            self.assertEqual(alignment.query.seq[43:94], "GLADTHTIEVTVDNEPVSLDITEESTSDLDKFNSGDKVTITYEKNDEGQLL")
            self.assertAlmostEqual(alignment.annotations["evalue"], 1e-04)
            self.assertAlmostEqual(alignment.annotations["bit score"], 31.6)
            self.assertAlmostEqual(alignment.score, 70)
            self.assertEqual(alignment.shape, (2, 59))
            self.assertAlmostEqual(alignment.annotations["% identity"], 33.90)
            self.assertEqual(alignment.annotations["identical"], 20)
            self.assertEqual(alignment.annotations["mismatches"], 31)
            self.assertEqual(alignment.annotations["positives"], 29)
            self.assertEqual(alignment.annotations["gap opens"], 1)
            self.assertEqual(alignment.annotations["gaps"], 8)
            self.assertAlmostEqual(alignment.annotations["% positives"], 49.15)
            self.assertEqual(alignment.annotations["query/sbjct frames"], "0/1")
            self.assertEqual(alignment.query.annotations["frame"], "0")
            self.assertEqual(alignment.target.annotations["frame"], "1")
            alignment = next(alignments)
            self.assertEqual(alignment.query.id, "gi|11464971:4-101")
            self.assertEqual(alignment.query.annotations["gi"], "0")
            self.assertEqual(alignment.query.annotations["acc."], "gi|11464971:4-101")
            self.assertEqual(alignment.query.annotations["acc.ver"], "gi|11464971:4-101")
            self.assertEqual(len(alignment.query.seq), 98)
            self.assertEqual(alignment.target.id, "gi|350596019|ref|XM_003360601.2|")
            self.assertEqual(alignment.target.annotations["ids"], "gi|350596019|ref|XM_003360601.2|")
            self.assertEqual(alignment.target.annotations["gi"], "0")
            self.assertEqual(alignment.target.annotations["gis"], "0")
            self.assertEqual(alignment.target.annotations["acc."], "gi|350596019|ref|XM_003360601.2|")
            self.assertEqual(alignment.target.annotations["acc.ver"], "gi|350596019|ref|XM_003360601.2|")
            self.assertEqual(alignment.target.annotations["length"], 772)
            self.assertTrue(
                numpy.array_equal(
                    alignment.coordinates, numpy.array([[0, 98], [0, 98]])
                )
            )
            self.assertEqual(alignment.target.annotations["start"], 94)
            self.assertEqual(alignment.target.annotations["end"], 388)
            self.assertEqual(alignment.target.seq, "KRIREGYLVKKGSMFNTWKPMWVILLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDGWVRDIKKAIK")
            self.assertEqual(alignment.query.seq[0:98], "KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK")
            self.assertAlmostEqual(alignment.annotations["evalue"], 2e-67)
            self.assertAlmostEqual(alignment.annotations["bit score"], 199)
            self.assertAlmostEqual(alignment.score, 506)
            self.assertEqual(alignment.shape, (2, 98))
            self.assertAlmostEqual(alignment.annotations["% identity"], 95.92)
            self.assertEqual(alignment.annotations["identical"], 94)
            self.assertEqual(alignment.annotations["mismatches"], 4)
            self.assertEqual(alignment.annotations["positives"], 96)
            self.assertEqual(alignment.annotations["gap opens"], 0)
            self.assertEqual(alignment.annotations["gaps"], 0)
            self.assertAlmostEqual(alignment.annotations["% positives"], 97.96)
            self.assertEqual(alignment.annotations["query/sbjct frames"], "0/2")
            self.assertEqual(alignment.query.annotations["frame"], "0")
            self.assertEqual(alignment.target.annotations["frame"], "2")
            alignment = next(alignments)
            self.assertEqual(alignment.query.id, "gi|11464971:4-101")
            self.assertEqual(alignment.query.annotations["gi"], "0")
            self.assertEqual(alignment.query.annotations["acc."], "gi|11464971:4-101")
            self.assertEqual(alignment.query.annotations["acc.ver"], "gi|11464971:4-101")
            self.assertEqual(len(alignment.query.seq), 98)
            self.assertEqual(alignment.target.id, "gi|350596019|ref|XM_003360601.2|")
            self.assertEqual(alignment.target.annotations["ids"], "gi|350596019|ref|XM_003360601.2|")
            self.assertEqual(alignment.target.annotations["gi"], "0")
            self.assertEqual(alignment.target.annotations["gis"], "0")
            self.assertEqual(alignment.target.annotations["acc."], "gi|350596019|ref|XM_003360601.2|")
            self.assertEqual(alignment.target.annotations["acc.ver"], "gi|350596019|ref|XM_003360601.2|")
            self.assertEqual(alignment.target.annotations["length"], 772)
            self.assertTrue(
                numpy.array_equal(
                    alignment.coordinates,
                    numpy.array([[ 0, 25, 26, 39, 42, 71],
                                 [29, 54, 54, 67, 67, 96]])
                )
            )
            self.assertEqual(alignment.target.annotations["start"], 541)
            self.assertEqual(alignment.target.annotations["end"], 754)
            self.assertEqual(alignment.target.seq, "LHYYDPAGGEDPLGAIHLRGCVVTSVESNTDGKNGFLWERAXXITADEVHYFLQAANPKERTEWIKAIQVA")
            self.assertEqual(alignment.query.seq[29:96], "IEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA")
            self.assertAlmostEqual(alignment.annotations["evalue"], 4e-05)
            self.assertAlmostEqual(alignment.annotations["bit score"], 32.7)
            self.assertAlmostEqual(alignment.score, 73)
            self.assertEqual(alignment.shape, (2, 71))
            self.assertAlmostEqual(alignment.annotations["% identity"], 29.58)
            self.assertEqual(alignment.annotations["identical"], 21)
            self.assertEqual(alignment.annotations["mismatches"], 46)
            self.assertEqual(alignment.annotations["positives"], 33)
            self.assertEqual(alignment.annotations["gap opens"], 2)
            self.assertEqual(alignment.annotations["gaps"], 4)
            self.assertAlmostEqual(alignment.annotations["% positives"], 46.48)
            self.assertEqual(alignment.annotations["query/sbjct frames"], "0/2")
            self.assertEqual(alignment.query.annotations["frame"], "0")
            self.assertEqual(alignment.target.annotations["frame"], "2")
            alignment = next(alignments)
            self.assertEqual(alignment.query.id, "gi|11464971:4-101")
            self.assertEqual(alignment.query.annotations["gi"], "0")
            self.assertEqual(alignment.query.annotations["acc."], "gi|11464971:4-101")
            self.assertEqual(alignment.query.annotations["acc.ver"], "gi|11464971:4-101")
            self.assertEqual(len(alignment.query.seq), 98)
            self.assertEqual(alignment.target.id, "gi|301779869|ref|XM_002925302.1|")
            self.assertEqual(alignment.target.annotations["ids"], "gi|301779869|ref|XM_002925302.1|")
            self.assertEqual(alignment.target.annotations["gi"], "0")
            self.assertEqual(alignment.target.annotations["gis"], "0")
            self.assertEqual(alignment.target.annotations["acc."], "gi|301779869|ref|XM_002925302.1|")
            self.assertEqual(alignment.target.annotations["acc.ver"], "gi|301779869|ref|XM_002925302.1|")
            self.assertEqual(alignment.target.annotations["length"], 1144)
            self.assertTrue(
                numpy.array_equal(
                    alignment.coordinates, numpy.array([[0, 98], [0, 98]]),
                )
            )
            self.assertEqual(alignment.target.annotations["start"], 77)
            self.assertEqual(alignment.target.annotations["end"], 371)
            self.assertEqual(alignment.target.seq, "KRIREGYLVKRGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK")
            self.assertEqual(alignment.query.seq[0:98], "KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK")
            self.assertAlmostEqual(alignment.annotations["evalue"], 2e-67)
            self.assertAlmostEqual(alignment.annotations["bit score"], 202)
            self.assertAlmostEqual(alignment.score, 515)
            self.assertEqual(alignment.shape, (2, 98))
            self.assertAlmostEqual(alignment.annotations["% identity"], 97.96)
            self.assertEqual(alignment.annotations["identical"], 96)
            self.assertEqual(alignment.annotations["mismatches"], 2)
            self.assertEqual(alignment.annotations["positives"], 97)
            self.assertEqual(alignment.annotations["gap opens"], 0)
            self.assertEqual(alignment.annotations["gaps"], 0)
            self.assertAlmostEqual(alignment.annotations["% positives"], 98.98)
            self.assertEqual(alignment.annotations["query/sbjct frames"], "0/3")
            self.assertEqual(alignment.query.annotations["frame"], "0")
            self.assertEqual(alignment.target.annotations["frame"], "3")
            alignment = next(alignments)
            self.assertEqual(alignment.query.id, "gi|11464971:4-101")
            self.assertEqual(alignment.query.annotations["gi"], "0")
            self.assertEqual(alignment.query.annotations["acc."], "gi|11464971:4-101")
            self.assertEqual(alignment.query.annotations["acc.ver"], "gi|11464971:4-101")
            self.assertEqual(len(alignment.query.seq), 98)
            self.assertEqual(alignment.target.id, "gi|301779869|ref|XM_002925302.1|")
            self.assertEqual(alignment.target.annotations["ids"], "gi|301779869|ref|XM_002925302.1|")
            self.assertEqual(alignment.target.annotations["gi"], "0")
            self.assertEqual(alignment.target.annotations["gis"], "0")
            self.assertEqual(alignment.target.annotations["acc."], "gi|301779869|ref|XM_002925302.1|")
            self.assertEqual(alignment.target.annotations["acc.ver"], "gi|301779869|ref|XM_002925302.1|")
            self.assertEqual(alignment.target.annotations["length"], 1144)
            self.assertTrue(
                numpy.array_equal(
                    alignment.coordinates,
                    numpy.array([[0, 27, 29, 54, 58, 100],
                                 [2, 29, 29, 54, 54,  96]])

                )
            )
            self.assertEqual(alignment.target.annotations["start"], 803)
            self.assertEqual(alignment.target.annotations["end"], 1103)
            self.assertEqual(alignment.target.seq, "IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGEDPLGAIHLRGCVVTSVESNPDVRKSEEENLFEIITADEVHYFLQAATPKERTEWIKAIQVA")
            self.assertEqual(alignment.query.seq[2:96], "IREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA")
            self.assertAlmostEqual(alignment.annotations["evalue"], 3e-09)
            self.assertAlmostEqual(alignment.annotations["bit score"], 45.1)
            self.assertAlmostEqual(alignment.score, 105)
            self.assertEqual(alignment.shape, (2, 100))
            self.assertAlmostEqual(alignment.annotations["% identity"], 30.00)
            self.assertEqual(alignment.annotations["identical"], 30)
            self.assertEqual(alignment.annotations["mismatches"], 64)
            self.assertEqual(alignment.annotations["positives"], 48)
            self.assertEqual(alignment.annotations["gap opens"], 2)
            self.assertEqual(alignment.annotations["gaps"], 6)
            self.assertAlmostEqual(alignment.annotations["% positives"], 48.00)
            self.assertEqual(alignment.annotations["query/sbjct frames"], "0/3")
            self.assertEqual(alignment.query.annotations["frame"], "0")
            self.assertEqual(alignment.target.annotations["frame"], "3")
            alignment = next(alignments)
            self.assertEqual(alignment.query.id, "gi|11464971:4-101")
            self.assertEqual(alignment.query.annotations["gi"], "0")
            self.assertEqual(alignment.query.annotations["acc."], "gi|11464971:4-101")
            self.assertEqual(alignment.query.annotations["acc.ver"], "gi|11464971:4-101")
            self.assertEqual(len(alignment.query.seq), 98)
            self.assertEqual(alignment.target.id, "gi|296223671|ref|XM_002757683.1|")
            self.assertEqual(alignment.target.annotations["ids"], "gi|296223671|ref|XM_002757683.1|")
            self.assertEqual(alignment.target.annotations["gi"], "0")
            self.assertEqual(alignment.target.annotations["gis"], "0")
            self.assertEqual(alignment.target.annotations["acc."], "gi|296223671|ref|XM_002757683.1|")
            self.assertEqual(alignment.target.annotations["acc.ver"], "gi|296223671|ref|XM_002757683.1|")
            self.assertEqual(alignment.target.annotations["length"], 1183)
            self.assertTrue(
                numpy.array_equal(
                    alignment.coordinates,
                    numpy.array([[0, 98], [0, 98]]),
                )
            )
            self.assertEqual(alignment.target.annotations["start"], 160)
            self.assertEqual(alignment.target.annotations["end"], 454)
            self.assertEqual(alignment.target.seq, "KRIREGYLVKKGSMFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK")
            self.assertEqual(alignment.query.seq[0:98], "KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK")
            self.assertAlmostEqual(alignment.annotations["evalue"], 4e-67)
            self.assertAlmostEqual(alignment.annotations["bit score"], 202)
            self.assertAlmostEqual(alignment.score, 515)
            self.assertEqual(alignment.shape, (2, 98))
            self.assertAlmostEqual(alignment.annotations["% identity"], 97.96)
            self.assertEqual(alignment.annotations["identical"], 96)
            self.assertEqual(alignment.annotations["mismatches"], 2)
            self.assertEqual(alignment.annotations["positives"], 97)
            self.assertEqual(alignment.annotations["gap opens"], 0)
            self.assertEqual(alignment.annotations["gaps"], 0)
            self.assertAlmostEqual(alignment.annotations["% positives"], 98.98)
            self.assertEqual(alignment.annotations["query/sbjct frames"], "0/2")
            self.assertEqual(alignment.query.annotations["frame"], "0")
            self.assertEqual(alignment.target.annotations["frame"], "2")
            alignment = next(alignments)
            self.assertEqual(alignment.query.id, "gi|11464971:4-101")
            self.assertEqual(alignment.query.annotations["gi"], "0")
            self.assertEqual(alignment.query.annotations["acc."], "gi|11464971:4-101")
            self.assertEqual(alignment.query.annotations["acc.ver"], "gi|11464971:4-101")
            self.assertEqual(len(alignment.query.seq), 98)
            self.assertEqual(alignment.target.id, "gi|296223671|ref|XM_002757683.1|")
            self.assertEqual(alignment.target.annotations["ids"], "gi|296223671|ref|XM_002757683.1|")
            self.assertEqual(alignment.target.annotations["gi"], "0")
            self.assertEqual(alignment.target.annotations["gis"], "0")
            self.assertEqual(alignment.target.annotations["acc."], "gi|296223671|ref|XM_002757683.1|")
            self.assertEqual(alignment.target.annotations["acc.ver"], "gi|296223671|ref|XM_002757683.1|")
            self.assertEqual(alignment.target.annotations["length"], 1183)
            self.assertTrue(
                numpy.array_equal(
                    alignment.coordinates,
                    numpy.array([[0, 27, 29, 64, 68, 100],
                                 [2, 29, 29, 64, 64,  96]])
                )
            )
            self.assertEqual(alignment.target.annotations["start"], 865)
            self.assertEqual(alignment.target.annotations["end"], 1165)
            self.assertEqual(alignment.target.seq, "IKQGCLLKQGHRRKNWKVRKFILREDPAYLHYYDPAGGEDPLGAIHLRGCVVTSVESNSDGRKSEEENLFEIITADEVHYFLQAATPKERTEWIKAIQVA")
            self.assertEqual(alignment.query.seq[2:96], "IREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA")
            self.assertAlmostEqual(alignment.annotations["evalue"], 3e-09)
            self.assertAlmostEqual(alignment.annotations["bit score"], 45.1)
            self.assertAlmostEqual(alignment.score, 105)
            self.assertEqual(alignment.shape, (2, 100))
            self.assertAlmostEqual(alignment.annotations["% identity"], 30.00)
            self.assertEqual(alignment.annotations["identical"], 30)
            self.assertEqual(alignment.annotations["mismatches"], 64)
            self.assertEqual(alignment.annotations["positives"], 48)
            self.assertEqual(alignment.annotations["gap opens"], 2)
            self.assertEqual(alignment.annotations["gaps"], 6)
            self.assertAlmostEqual(alignment.annotations["% positives"], 48.00)
            self.assertEqual(alignment.annotations["query/sbjct frames"], "0/2")
            self.assertEqual(alignment.query.annotations["frame"], "0")
            self.assertEqual(alignment.target.annotations["frame"], "2")
            alignment = next(alignments)
            self.assertEqual(alignment.query.id, "gi|11464971:4-101")
            self.assertEqual(alignment.query.annotations["gi"], "0")
            self.assertEqual(alignment.query.annotations["acc."], "gi|11464971:4-101")
            self.assertEqual(alignment.query.annotations["acc.ver"], "gi|11464971:4-101")
            self.assertEqual(len(alignment.query.seq), 98)
            self.assertEqual(alignment.target.id, "gi|338714227|ref|XM_001492113.3|")
            self.assertEqual(alignment.target.annotations["ids"], "gi|338714227|ref|XM_001492113.3|")
            self.assertEqual(alignment.target.annotations["gi"], "0")
            self.assertEqual(alignment.target.annotations["gis"], "0")
            self.assertEqual(alignment.target.annotations["acc."], "gi|338714227|ref|XM_001492113.3|")
            self.assertEqual(alignment.target.annotations["acc.ver"], "gi|338714227|ref|XM_001492113.3|")
            self.assertEqual(alignment.target.annotations["length"], 1390)
            self.assertTrue(
                numpy.array_equal(
                    alignment.coordinates,
                    numpy.array([[0, 98], [0, 98]]),
                )
            )
            self.assertEqual(alignment.target.annotations["start"], 172)
            self.assertEqual(alignment.target.annotations["end"], 466)
            self.assertEqual(alignment.target.seq, "KRIREGYLVKRGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVFKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK")
            self.assertEqual(alignment.query.seq[0:98], "KRIREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKAIK")
            self.assertAlmostEqual(alignment.annotations["evalue"], 2e-66)
            self.assertAlmostEqual(alignment.annotations["bit score"], 202)
            self.assertAlmostEqual(alignment.score, 515)
            self.assertEqual(alignment.shape, (2, 98))
            self.assertAlmostEqual(alignment.annotations["% identity"], 97.96)
            self.assertEqual(alignment.annotations["identical"], 96)
            self.assertEqual(alignment.annotations["mismatches"], 2)
            self.assertEqual(alignment.annotations["positives"], 97)
            self.assertEqual(alignment.annotations["gap opens"], 0)
            self.assertEqual(alignment.annotations["gaps"], 0)
            self.assertAlmostEqual(alignment.annotations["% positives"], 98.98)
            self.assertEqual(alignment.annotations["query/sbjct frames"], "0/2")
            self.assertEqual(alignment.query.annotations["frame"], "0")
            self.assertEqual(alignment.target.annotations["frame"], "2")
            alignment = next(alignments)
            self.assertEqual(alignment.query.id, "gi|11464971:4-101")
            self.assertEqual(alignment.query.annotations["gi"], "0")
            self.assertEqual(alignment.query.annotations["acc."], "gi|11464971:4-101")
            self.assertEqual(alignment.query.annotations["acc.ver"], "gi|11464971:4-101")
            self.assertEqual(len(alignment.query.seq), 98)
            self.assertEqual(alignment.target.id, "gi|338714227|ref|XM_001492113.3|")
            self.assertEqual(alignment.target.annotations["ids"], "gi|338714227|ref|XM_001492113.3|")
            self.assertEqual(alignment.target.annotations["gi"], "0")
            self.assertEqual(alignment.target.annotations["gis"], "0")
            self.assertEqual(alignment.target.annotations["acc."], "gi|338714227|ref|XM_001492113.3|")
            self.assertEqual(alignment.target.annotations["acc.ver"], "gi|338714227|ref|XM_001492113.3|")
            self.assertEqual(alignment.target.annotations["length"], 1390)
            self.assertTrue(
                numpy.array_equal(
                    alignment.coordinates,
                    numpy.array([[0, 27, 29, 54, 58, 100],
                                 [2, 29, 29, 54, 54,  96]])
                )
            )
            self.assertEqual(alignment.target.annotations["start"], 898)
            self.assertEqual(alignment.target.annotations["end"], 1198)
            self.assertEqual(alignment.target.seq, "IKQGCLLKQGHRRKNWKVRKFVLREDPAYVHYYDPAGGEEPLGAIHLRGCVVTSVEGNPDGKKSEEENLFEIITADEVHYFLQAATPKERTEWIKAIQVA")
            self.assertEqual(alignment.query.seq[2:96], "IREGYLVKKGSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTSPCQDFGKRMFVLKITTTKQQDHFFQAAFLEERDAWVRDIKKA")
            self.assertAlmostEqual(alignment.annotations["evalue"], 1e-09)
            self.assertAlmostEqual(alignment.annotations["bit score"], 46.6)
            self.assertAlmostEqual(alignment.score, 109)
            self.assertEqual(alignment.shape, (2, 100))
            self.assertAlmostEqual(alignment.annotations["% identity"], 31.00)
            self.assertEqual(alignment.annotations["identical"], 31)
            self.assertEqual(alignment.annotations["mismatches"], 63)
            self.assertEqual(alignment.annotations["positives"], 48)
            self.assertEqual(alignment.annotations["gap opens"], 2)
            self.assertEqual(alignment.annotations["gaps"], 6)
            self.assertAlmostEqual(alignment.annotations["% positives"], 48.00)
            self.assertEqual(alignment.annotations["query/sbjct frames"], "0/2")
            self.assertEqual(alignment.query.annotations["frame"], "0")
            self.assertEqual(alignment.target.annotations["frame"], "2")
            alignment = next(alignments)
            self.assertEqual(alignment.query.id, "gi|11464971:4-101")
            self.assertEqual(alignment.query.annotations["gi"], "0")
            self.assertEqual(alignment.query.annotations["acc."], "gi|11464971:4-101")
            self.assertEqual(alignment.query.annotations["acc.ver"], "gi|11464971:4-101")
            self.assertEqual(len(alignment.query.seq), 98)
            self.assertEqual(alignment.target.id, "gi|365982352|ref|XM_003667962.1|")
            self.assertEqual(alignment.target.annotations["ids"], "gi|365982352|ref|XM_003667962.1|")
            self.assertEqual(alignment.target.annotations["gi"], "0")
            self.assertEqual(alignment.target.annotations["gis"], "0")
            self.assertEqual(alignment.target.annotations["acc."], "gi|365982352|ref|XM_003667962.1|")
            self.assertEqual(alignment.target.annotations["acc.ver"], "gi|365982352|ref|XM_003667962.1|")
            self.assertEqual(alignment.target.annotations["length"], 4932)
            self.assertTrue(
                numpy.array_equal(
                    alignment.coordinates,
                    numpy.array([[ 0, 15, 24, 52],
                                 [11, 26, 26, 54]])
                )
            )
            self.assertEqual(alignment.target.annotations["start"], 3180)
            self.assertEqual(alignment.target.annotations["end"], 3336)
            self.assertEqual(alignment.target.seq, "GSCFPTWDLIFIEVLNPFLKEKLWEADNEEISKFVDLTLKGLVDLYPSHFTS")
            self.assertEqual(alignment.query.seq[11:54], "GSVFNTWKPMWVVLLEDGIEFYKKKSDNSPKGMIPLKGSTLTS")
            self.assertAlmostEqual(alignment.annotations["evalue"], 1.7)
            self.assertAlmostEqual(alignment.annotations["bit score"], 19.6)
            self.assertAlmostEqual(alignment.score, 39)
            self.assertEqual(alignment.shape, (2, 52))
            self.assertAlmostEqual(alignment.annotations["% identity"], 30.77)
            self.assertEqual(alignment.annotations["identical"], 16)
            self.assertEqual(alignment.annotations["mismatches"], 27)
            self.assertEqual(alignment.annotations["positives"], 23)
            self.assertEqual(alignment.annotations["gap opens"], 1)
            self.assertEqual(alignment.annotations["gaps"], 9)
            self.assertAlmostEqual(alignment.annotations["% positives"], 44.23)
            self.assertEqual(alignment.annotations["query/sbjct frames"], "0/1")
            self.assertEqual(alignment.query.annotations["frame"], "0")
            self.assertEqual(alignment.target.annotations["frame"], "1")
            with self.assertRaises(StopIteration):
                next(alignments)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
