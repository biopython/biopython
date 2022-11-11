# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for Bio.Align Exonerate parsers."""

import os
import io
import unittest

from Bio import Align


try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.exonerate."
    ) from None


class Exonerate_est2genome(unittest.TestCase):
    def test_exn_22_m_est2genome_cigar(self):
        """Test parsing exn_22_m_est2genome_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_est2genome_cigar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)

    def check_cigar(self, alignments):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m est2genome ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar yes --showvulgar no",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 6150)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 0 1230 + gi|330443520|ref|NC_001136.10| 1319275 1318045 - 6150 M 1230
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 439)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 85010,  85021,  85021,  85036,  85036,  85040,
                               85040,  85041,  85041,  85049,  85049,  85066,
                              253974, 253978, 253979, 253987, 253987, 253990,
                              253990, 254023, 254024, 254031, 254033, 254135,
                              350959, 350973, 350975, 350985, 350985, 350990,
                              350992, 351002, 351002, 351006, 351007, 351027,
                              351027, 351042, 351043, 351048, 351048, 351052,
                              473170, 473190, 473195, 473201],
                             [     0,     11,     12,     27,     29,     33,
                                  34,     35,     36,     44,     48,     65,
                                  65,     69,     69,     77,     78,     81,
                                  83,    116,    116,    123,    123,    225,
                                 225,    239,    239,    249,    251,    256,
                                 256,    266,    268,    272,    272,    292,
                                 293,    308,    308,    313,    316,    320,
                                 320,    340,    340,    346]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 0 346 + gi|330443688|ref|NC_001145.3| 85010 473201 + 439 M 11 I 1 M 15 I 2 M 4 I 1 M 1 I 1 M 8 I 4 M 17 D 168908 M 4 D 1 M 8 I 1 M 3 I 2 M 33 D 1 M 7 D 2 M 102 D 96824 M 14 D 2 M 10 I 2 M 5 D 2 M 10 I 2 M 4 D 1 M 20 I 1 M 15 D 1 M 5 I 3 M 4 D 122118 M 20 D 5 M 6
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 263)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[130198, 130184, 130183, 130179, 130179, 130154,
                              130153, 130144, 130138, 130096, 130096, 130080,
                              130078, 130071, 130070, 130067, 130067, 130044,
                              130044, 130038, 120681, 120680, 120680, 120669,
                              120668, 120656, 120656, 120647, 120646, 120636,
                              120636, 120618, 120617, 120612,  11487,  11471,
                               11471,  11467,  11467,  11456,  11456,  11448,
                               11448,  11426,  11424,  11420,  11418,  11384,
                               11383,  11380,  11379,  11338],
                             [    25,     39,     39,     43,     45,     70,
                                  70,     79,     79,    121,    123,    139,
                                 139,    146,    146,    149,    151,    174,
                                 177,    183,    183,    184,    185,    196,
                                 196,    208,    209,    218,    218,    228,
                                 229,    247,    247,    252,    252,    268,
                                 272,    276,    277,    288,    293,    301,
                                 302,    324,    324,    328,    328,    362,
                                 362,    365,    365,    406]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 25 406 + gi|330443688|ref|NC_001145.3| 130198 11338 - 263 M 14 D 1 M 4 I 2 M 25 D 1 M 9 D 6 M 42 I 2 M 16 D 2 M 7 D 1 M 3 I 2 M 23 I 3 M 6 D 9357 M 1 I 1 M 11 D 1 M 12 I 1 M 9 D 1 M 10 I 1 M 18 D 1 M 5 D 109125 M 16 I 4 M 4 I 1 M 11 I 5 M 8 I 1 M 22 D 2 M 4 D 2 M 34 D 1 M 3 D 1 M 41
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_est2genome_vulgar(self):
        """Test parsing exn_22_m_est2genome_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_est2genome_vulgar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments, check_operations=False)

    def check_vulgar(self, alignments, check_operations=True):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m est2genome ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar no --showvulgar yes",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 6150)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"M"))
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: gi|296143771|ref|NM_001180731.1| 0 1230 + gi|330443520|ref|NC_001136.10| 1319275 1318045 - 6150 M 1230 1230
""",
            )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 439)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 85010,  85021,  85021,  85036,  85036,  85040,
                               85040,  85041,  85041,  85049,  85049,  85066,
                               85068, 253972, 253974, 253978, 253979, 253987,
                              253987, 253990, 253990, 254023, 254024, 254031,
                              254033, 254135, 254137, 350957, 350959, 350973,
                              350975, 350985, 350985, 350990, 350992, 351002,
                              351002, 351006, 351007, 351027, 351027, 351042,
                              351043, 351048, 351048, 351052, 351054, 473168,
                              473170, 473190, 473195, 473201],
                             [     0,     11,     12,     27,     29,     33,
                                  34,     35,     36,     44,     48,     65,
                                  65,     65,     65,     69,     69,     77,
                                  78,     81,     83,    116,    116,    123,
                                 123,    225,    225,    225,    225,    239,
                                 239,    249,    251,    256,    256,    266,
                                 268,    272,    272,    292,    293,    308,
                                 308,    313,    316,    320,    320,    320,
                                 320,    340,    340,    346]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(
                alignment.operations,
                bytearray(b"MIMIMIMIMIM5N3MDMIMIMDMDM5N3MDMIMDMIMDMIMDMIM5N3MDM"),
            )
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: gi|296143771|ref|NM_001180731.1| 0 346 + gi|330443688|ref|NC_001145.3| 85010 473201 + 439 M 11 11 G 1 0 M 15 15 G 2 0 M 4 4 G 1 0 M 1 1 G 1 0 M 8 8 G 4 0 M 17 17 5 0 2 I 0 168904 3 0 2 M 4 4 G 0 1 M 8 8 G 1 0 M 3 3 G 2 0 M 33 33 G 0 1 M 7 7 G 0 2 M 102 102 5 0 2 I 0 96820 3 0 2 M 14 14 G 0 2 M 10 10 G 2 0 M 5 5 G 0 2 M 10 10 G 2 0 M 4 4 G 0 1 M 20 20 G 1 0 M 15 15 G 0 1 M 5 5 G 3 0 M 4 4 5 0 2 I 0 122114 3 0 2 M 20 20 G 0 5 M 6 6
""",
            )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 263)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[130198, 130184, 130183, 130179, 130179, 130154,
                              130153, 130144, 130138, 130096, 130096, 130080,
                              130078, 130071, 130070, 130067, 130067, 130044,
                              130044, 130038, 130036, 120683, 120681, 120680,
                              120680, 120669, 120668, 120656, 120656, 120647,
                              120646, 120636, 120636, 120618, 120617, 120612,
                              120610,  11489,  11487,  11471,  11471,  11467,
                               11467,  11456,  11456,  11448,  11448,  11426,
                               11424,  11420,  11418,  11384,  11383,  11380,
                               11379,  11338],
                             [    25,     39,     39,     43,     45,     70,
                                  70,     79,     79,    121,    123,    139,
                                 139,    146,    146,    149,    151,    174,
                                 177,    183,    183,    183,    183,    184,
                                 185,    196,    196,    208,    209,    218,
                                 218,    228,    229,    247,    247,    252,
                                 252,    252,    252,    268,    272,    276,
                                 277,    288,    293,    301,    302,    324,
                                 324,    328,    328,    362,    362,    365,
                                 365,    406]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(
                alignment.operations,
                bytearray(b"MDMIMDMDMIMDMDMIMIM3N5MIMDMIMDMIMDM3N5MIMIMIMIMDMDMDMDM"),
            )
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: gi|296143771|ref|NM_001180731.1| 25 406 + gi|330443688|ref|NC_001145.3| 130198 11338 - 263 M 14 14 G 0 1 M 4 4 G 2 0 M 25 25 G 0 1 M 9 9 G 0 6 M 42 42 G 2 0 M 16 16 G 0 2 M 7 7 G 0 1 M 3 3 G 2 0 M 23 23 G 3 0 M 6 6 3 0 2 I 0 9353 5 0 2 M 1 1 G 1 0 M 11 11 G 0 1 M 12 12 G 1 0 M 9 9 G 0 1 M 10 10 G 1 0 M 18 18 G 0 1 M 5 5 3 0 2 I 0 109121 5 0 2 M 16 16 G 4 0 M 4 4 G 1 0 M 11 11 G 5 0 M 8 8 G 1 0 M 22 22 G 0 2 M 4 4 G 0 2 M 34 34 G 0 1 M 3 3 G 0 1 M 41 41
""",
            )
        with self.assertRaises(StopIteration):
            next(alignments)


class Exonerate_affine_local(unittest.TestCase):
    def test_exn_22_m_affine_local_cigar(self):
        """Test parsing exn_22_m_affine_local_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_affine_local_cigar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)

    def check_cigar(self, alignments):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m affine:local ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar yes --showvulgar no",
        )
        self.assertEqual(alignments.metadata["Hostname"], "Michiels-MacBook-Pro.local")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 6150)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 0 1230 + gi|330443520|ref|NC_001136.10| 1319275 1318045 - 6150 M 1230
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 359)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[253990, 254024, 254025, 254031, 254033, 254134,
                              254136, 254142, 254142, 254145, 254145, 254157,
                              254158, 254167, 254167, 254193, 254194, 254228,
                              254233, 254246, 254248, 254255, 254256, 254266,
                              254267, 254281, 254284, 254293, 254294, 254314,
                              254316, 254318, 254319, 254326, 254327, 254329,
                              254330, 254336, 254337, 254340, 254340, 254342,
                              254342, 254352, 254353, 254355, 254358, 254367,
                              254367, 254376, 254376, 254383, 254384, 254390,
                              254390, 254398, 254399, 254431, 254431, 254440,
                              254440, 254449, 254451, 254474],
                             [    83,    117,    117,    123,    123,    224,
                                 224,    230,    232,    235,    236,    248,
                                 248,    257,    258,    284,    284,    318,
                                 318,    331,    331,    338,    338,    348,
                                 348,    362,    362,    371,    371,    391,
                                 391,    393,    393,    400,    400,    402,
                                 402,    408,    408,    411,    412,    414,
                                 417,    427,    427,    429,    429,    438,
                                 440,    449,    452,    459,    459,    465,
                                 467,    475,    475,    507,    509,    518,
                                 520,    529,    529,    552]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 83 552 + gi|330443688|ref|NC_001145.3| 253990 254474 + 359 M 34 D 1 M 6 D 2 M 101 D 2 M 6 I 2 M 3 I 1 M 12 D 1 M 9 I 1 M 26 D 1 M 34 D 5 M 13 D 2 M 7 D 1 M 10 D 1 M 14 D 3 M 9 D 1 M 20 D 2 M 2 D 1 M 7 D 1 M 2 D 1 M 6 D 1 M 3 I 1 M 2 I 3 M 10 D 1 M 2 D 3 M 9 I 2 M 9 I 3 M 7 D 1 M 6 I 2 M 8 D 1 M 32 I 2 M 9 I 2 M 9 D 2 M 23
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443715|ref|NC_001146.8|")
        self.assertEqual(alignment.score, 219)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[454073, 454087, 454087, 454094, 454094, 454097,
                              454097, 454109, 454110, 454133, 454134, 454150,
                              454150, 454167, 454169, 454175, 454176, 454183,
                              454184, 454192, 454192, 454204, 454205, 454208,
                              454209, 454223, 454224, 454234, 454235, 454255,
                              454256, 454275, 454277, 454286, 454286, 454294,
                              454294, 454335, 454335, 454377, 454378, 454386,
                              454389, 454390, 454391, 454417, 454417, 454428,
                              454428, 454440, 454441, 454459, 454460, 454478,
                              454479, 454493, 454494, 454524, 454524, 454531],
                             [    60,     74,     75,     82,     83,     86,
                                  87,     99,     99,    122,    122,    138,
                                 140,    157,    157,    163,    163,    170,
                                 170,    178,    179,    191,    191,    194,
                                 194,    208,    208,    218,    218,    238,
                                 238,    257,    257,    266,    271,    279,
                                 280,    321,    324,    366,    366,    374,
                                 374,    375,    375,    401,    403,    414,
                                 417,    429,    429,    447,    447,    465,
                                 465,    479,    479,    509,    510,    517]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 60 517 + gi|330443715|ref|NC_001146.8| 454073 454531 + 219 M 14 I 1 M 7 I 1 M 3 I 1 M 12 D 1 M 23 D 1 M 16 I 2 M 17 D 2 M 6 D 1 M 7 D 1 M 8 I 1 M 12 D 1 M 3 D 1 M 14 D 1 M 10 D 1 M 20 D 1 M 19 D 2 M 9 I 5 M 8 I 1 M 41 I 3 M 42 D 1 M 8 D 3 M 1 D 1 M 26 I 2 M 11 I 3 M 12 D 1 M 18 D 1 M 18 D 1 M 14 D 1 M 30 I 1 M 7
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_affine_local_vulgar(self):
        """Test parsing exn_22_m_affine_local_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_affine_local_vulgar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments, check_operations=False)

    def check_vulgar(self, alignments, check_operations=True):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m affine:local ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar no --showvulgar yes",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 6150)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"M"))
        self.assertEqual(
            alignment.format("exonerate", "vulgar"),
            """\
vulgar: gi|296143771|ref|NM_001180731.1| 0 1230 + gi|330443520|ref|NC_001136.10| 1319275 1318045 - 6150 M 1230 1230
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 359)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[253990, 254024, 254025, 254031, 254033, 254134,
                              254136, 254142, 254142, 254145, 254145, 254157,
                              254158, 254167, 254167, 254193, 254194, 254228,
                              254233, 254246, 254248, 254255, 254256, 254266,
                              254267, 254281, 254284, 254293, 254294, 254314,
                              254316, 254318, 254319, 254326, 254327, 254329,
                              254330, 254336, 254337, 254340, 254340, 254342,
                              254342, 254352, 254353, 254355, 254358, 254367,
                              254367, 254376, 254376, 254383, 254384, 254390,
                              254390, 254398, 254399, 254431, 254431, 254440,
                              254440, 254449, 254451, 254474],
                             [    83,    117,    117,    123,    123,    224,
                                 224,    230,    232,    235,    236,    248,
                                 248,    257,    258,    284,    284,    318,
                                 318,    331,    331,    338,    338,    348,
                                 348,    362,    362,    371,    371,    391,
                                 391,    393,    393,    400,    400,    402,
                                 402,    408,    408,    411,    412,    414,
                                 417,    427,    427,    429,    429,    438,
                                 440,    449,    452,    459,    459,    465,
                                 467,    475,    475,    507,    509,    518,
                                 520,    529,    529,    552]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(
                alignment.operations,
                bytearray(
                    b"MDMDMDMIMIMDMIMDMDMDMDMDMDMDMDMDMDMDMDMIMIMDMDMIMIMDMIMDMIMIMDM"
                ),
            )
        self.assertEqual(
            alignment.format("exonerate", "vulgar"),
            """\
vulgar: gi|296143771|ref|NM_001180731.1| 83 552 + gi|330443688|ref|NC_001145.3| 253990 254474 + 359 M 34 34 G 0 1 M 6 6 G 0 2 M 101 101 G 0 2 M 6 6 G 2 0 M 3 3 G 1 0 M 12 12 G 0 1 M 9 9 G 1 0 M 26 26 G 0 1 M 34 34 G 0 5 M 13 13 G 0 2 M 7 7 G 0 1 M 10 10 G 0 1 M 14 14 G 0 3 M 9 9 G 0 1 M 20 20 G 0 2 M 2 2 G 0 1 M 7 7 G 0 1 M 2 2 G 0 1 M 6 6 G 0 1 M 3 3 G 1 0 M 2 2 G 3 0 M 10 10 G 0 1 M 2 2 G 0 3 M 9 9 G 2 0 M 9 9 G 3 0 M 7 7 G 0 1 M 6 6 G 2 0 M 8 8 G 0 1 M 32 32 G 2 0 M 9 9 G 2 0 M 9 9 G 0 2 M 23 23
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443715|ref|NC_001146.8|")
        self.assertEqual(alignment.score, 219)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[454073, 454087, 454087, 454094, 454094, 454097,
                              454097, 454109, 454110, 454133, 454134, 454150,
                              454150, 454167, 454169, 454175, 454176, 454183,
                              454184, 454192, 454192, 454204, 454205, 454208,
                              454209, 454223, 454224, 454234, 454235, 454255,
                              454256, 454275, 454277, 454286, 454286, 454294,
                              454294, 454335, 454335, 454377, 454378, 454386,
                              454389, 454390, 454391, 454417, 454417, 454428,
                              454428, 454440, 454441, 454459, 454460, 454478,
                              454479, 454493, 454494, 454524, 454524, 454531],
                             [    60,     74,     75,     82,     83,     86,
                                  87,     99,     99,    122,    122,    138,
                                 140,    157,    157,    163,    163,    170,
                                 170,    178,    179,    191,    191,    194,
                                 194,    208,    208,    218,    218,    238,
                                 238,    257,    257,    266,    271,    279,
                                 280,    321,    324,    366,    366,    374,
                                 374,    375,    375,    401,    403,    414,
                                 417,    429,    429,    447,    447,    465,
                                 465,    479,    479,    509,    510,    517]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(
                alignment.operations,
                bytearray(
                    b"MIMIMIMDMDMIMDMDMDMIMDMDMDMDMDMDMIMIMIMDMDMDMIMIMDMDMDMDMIM"
                ),
            )
        self.assertEqual(
            alignment.format("exonerate", "vulgar"),
            """\
vulgar: gi|296143771|ref|NM_001180731.1| 60 517 + gi|330443715|ref|NC_001146.8| 454073 454531 + 219 M 14 14 G 1 0 M 7 7 G 1 0 M 3 3 G 1 0 M 12 12 G 0 1 M 23 23 G 0 1 M 16 16 G 2 0 M 17 17 G 0 2 M 6 6 G 0 1 M 7 7 G 0 1 M 8 8 G 1 0 M 12 12 G 0 1 M 3 3 G 0 1 M 14 14 G 0 1 M 10 10 G 0 1 M 20 20 G 0 1 M 19 19 G 0 2 M 9 9 G 5 0 M 8 8 G 1 0 M 41 41 G 3 0 M 42 42 G 0 1 M 8 8 G 0 3 M 1 1 G 0 1 M 26 26 G 2 0 M 11 11 G 3 0 M 12 12 G 0 1 M 18 18 G 0 1 M 18 18 G 0 1 M 14 14 G 0 1 M 30 30 G 1 0 M 7 7
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)


class Exonerate_cdna2genome(unittest.TestCase):
    def test_exn_22_m_cdna2genome_cigar(self):
        """Test parsing exn_22_m_cdna2genome_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_cdna2genome_cigar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)

    def check_cigar(self, alignments):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m cdna2genome ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar yes --showvulgar no",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 6146)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                # fmt: off
# flake8: noqa
                alignment.coordinates,
                numpy.array([[1319275, 1319274, 1319271, 1318045],
                             [      0,       1,       4,    1230]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 0 1230 + gi|330443520|ref|NC_001136.10| 1319275 1318045 - 6146 M 1 M 3 M 1226
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 6146)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[1318045, 1318174, 1318177, 1319275],
                             [   1230,    1101,    1098,       0]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 1230 0 - gi|330443520|ref|NC_001136.10| 1318045 1319275 + 6146 M 129 M 3 M 1098
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 518)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 85010,  85021,  85021,  85036,  85036,  85040,
                               85040,  85041,  85041,  85049,  85049,  85066,
                              253974, 253978, 253979, 253987, 253987, 253990,
                              253990, 254023, 254025, 254032, 254033, 254135,
                              350959, 350973, 350975, 350985, 350985, 350990,
                              350992, 351002, 351002, 351006, 351007, 351027,
                              351027, 351042, 351043, 351048, 351048, 351052,
                              473170, 473190, 473195, 473201, 667040, 667052,
                              667054, 667059, 667059, 667066, 667068, 667069,
                              667070, 667082, 667157, 667163, 667163, 667167,
                              667168, 667170, 667171, 667174, 667175, 667216],
                             [     0,     11,     12,     27,     29,     33,
                                  34,     35,     36,     44,     48,     65,
                                  65,     69,     69,     77,     79,     82,
                                  83,    116,    116,    123,    123,    225,
                                 225,    239,    239,    249,    251,    256,
                                 256,    266,    268,    272,    272,    292,
                                 293,    308,    308,    313,    316,    320,
                                 320,    340,    340,    346,    346,    358,
                                 358,    363,    364,    371,    371,    372,
                                 372,    384,    459,    465,    466,    470,
                                 470,    472,    472,    475,    475,    516]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 0 516 + gi|330443688|ref|NC_001145.3| 85010 667216 + 518 M 11 I 1 M 15 I 2 M 4 I 1 M 1 I 1 M 8 I 4 M 17 D 168908 M 4 D 1 M 8 I 2 M 3 I 1 M 33 D 2 M 7 D 1 M 102 D 96824 M 14 D 2 M 10 I 2 M 5 D 2 M 10 I 2 M 4 D 1 M 20 I 1 M 15 D 1 M 5 I 3 M 4 D 122118 M 20 D 5 M 6 D 193839 M 12 D 2 M 5 I 1 M 7 D 2 M 1 D 1 M 12 M 75 M 6 I 1 M 4 D 1 M 2 D 1 M 3 D 1 M 41
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_cdna2genome_vulgar(self):
        """Test parsing exn_22_m_cdna2genome_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_cdna2genome_vulgar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments, check_operations=False)

    def check_vulgar(self, alignments, check_operations=True):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m cdna2genome ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar no --showvulgar yes",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 6146)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[1319275, 1319274, 1319271, 1318045],
                             [      0,       1,       4,    1230]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"MCM"))
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: gi|296143771|ref|NM_001180731.1| 0 1230 + gi|330443520|ref|NC_001136.10| 1319275 1318045 - 6146 M 1 1 C 3 3 M 1226 1226
""",
            )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 6146)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[1318045, 1318174, 1318177, 1319275],
                             [   1230,    1101,    1098,       0]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"MCM"))
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: gi|296143771|ref|NM_001180731.1| 1230 0 - gi|330443520|ref|NC_001136.10| 1318045 1319275 + 6146 M 129 129 C 3 3 M 1098 1098
""",
            )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 518)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 85010,  85021,  85021,  85036,  85036,  85040,
                               85040,  85041,  85041,  85049,  85049,  85066,
                               85068, 253972, 253974, 253978, 253979, 253987,
                              253987, 253990, 253990, 254023, 254025, 254032,
                              254033, 254135, 254137, 350957, 350959, 350973,
                              350975, 350985, 350985, 350990, 350992, 351002,
                              351002, 351006, 351007, 351027, 351027, 351042,
                              351043, 351048, 351048, 351052, 351054, 473168,
                              473170, 473190, 473195, 473201, 473203, 667038,
                              667040, 667052, 667054, 667059, 667059, 667066,
                              667068, 667069, 667070, 667082, 667157, 667163,
                              667163, 667167, 667168, 667170, 667171, 667174,
                              667175, 667216],
                             [     0,     11,     12,     27,     29,     33,
                                  34,     35,     36,     44,     48,     65,
                                  65,     65,     65,     69,     69,     77,
                                  79,     82,     83,    116,    116,    123,
                                 123,    225,    225,    225,    225,    239,
                                 239,    249,    251,    256,    256,    266,
                                 268,    272,    272,    292,    293,    308,
                                 308,    313,    316,    320,    320,    320,
                                 320,    340,    340,    346,    346,    346,
                                 346,    358,    358,    363,    364,    371,
                                 371,    372,    372,    384,    459,    465,
                                 466,    470,    470,    472,    472,    475,
                                 475,    516]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(
                alignment.operations,
                bytearray(
                    b"MIMIMIMIMIM5N3MDMIMIMDMDM5N3MDMIMDMIMDMIMDMIM5N3MDM5N3MDMIMDMDMCMIMDMDMDM"
                ),
            )
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: gi|296143771|ref|NM_001180731.1| 0 516 + gi|330443688|ref|NC_001145.3| 85010 667216 + 518 M 11 11 G 1 0 M 15 15 G 2 0 M 4 4 G 1 0 M 1 1 G 1 0 M 8 8 G 4 0 M 17 17 5 0 2 I 0 168904 3 0 2 M 4 4 G 0 1 M 8 8 G 2 0 M 3 3 G 1 0 M 33 33 G 0 2 M 7 7 G 0 1 M 102 102 5 0 2 I 0 96820 3 0 2 M 14 14 G 0 2 M 10 10 G 2 0 M 5 5 G 0 2 M 10 10 G 2 0 M 4 4 G 0 1 M 20 20 G 1 0 M 15 15 G 0 1 M 5 5 G 3 0 M 4 4 5 0 2 I 0 122114 3 0 2 M 20 20 G 0 5 M 6 6 5 0 2 I 0 193835 3 0 2 M 12 12 G 0 2 M 5 5 G 1 0 M 7 7 G 0 2 M 1 1 G 0 1 M 12 12 C 75 75 M 6 6 G 1 0 M 4 4 G 0 1 M 2 2 G 0 1 M 3 3 G 0 1 M 41 41
""",
            )
        with self.assertRaises(StopIteration):
            next(alignments)


class Exonerate_coding2coding(unittest.TestCase):
    def test_exn_22_m_coding2coding_cigar(self):
        """Test parsing exn_22_m_coding2coding_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_coding2coding_cigar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)

    def check_cigar(self, alignments):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m coding2coding ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar yes --showvulgar no",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2151)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1318047, 1319274], [1228, 1]])
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 1228 1 - gi|330443520|ref|NC_001136.10| 1318047 1319274 + 2151 M 1227
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2106)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 0 1230 + gi|330443520|ref|NC_001136.10| 1319275 1318045 - 2106 M 1230
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 116)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[255638, 255743, 255743, 255794],
                             [  1065,   1170,   1173,   1224]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 1065 1224 + gi|330443688|ref|NC_001145.3| 255638 255794 + 116 M 105 I 3 M 51
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_coding2coding_vulgar(self):
        """Test parsing exn_22_m_coding2coding_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_coding2coding_vulgar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments, check_operations=False)

    def check_vulgar(self, alignments, check_operations=True):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m coding2coding ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar no --showvulgar yes",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2151)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1318047, 1319274], [1228, 1]])
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"C"))
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: gi|296143771|ref|NM_001180731.1| 1228 1 - gi|330443520|ref|NC_001136.10| 1318047 1319274 + 2151 C 1227 1227
""",
            )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2106)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"C"))
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: gi|296143771|ref|NM_001180731.1| 0 1230 + gi|330443520|ref|NC_001136.10| 1319275 1318045 - 2106 C 1230 1230
""",
            )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 116)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[255638, 255743, 255743, 255794],
                             [  1065,   1170,   1173,   1224]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"CIC"))
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: gi|296143771|ref|NM_001180731.1| 1065 1224 + gi|330443688|ref|NC_001145.3| 255638 255794 + 116 C 105 105 G 3 0 C 51 51
""",
            )
        with self.assertRaises(StopIteration):
            next(alignments)


class Exonerate_coding2genome(unittest.TestCase):
    def test_exn_22_m_coding2genome_cigar(self):
        """Test parsing exn_22_m_coding2genome_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_coding2genome_cigar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)

    def check_cigar(self, alignments):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m coding2genome ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar yes --showvulgar no",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2151)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                # fmt: off
# flake8: noqa
                alignment.coordinates, numpy.array([[1318047, 1319274], [1228, 1]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 1228 1 - gi|330443520|ref|NC_001136.10| 1318047 1319274 + 2151 M 1227
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2106)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 0 1230 + gi|330443520|ref|NC_001136.10| 1319275 1318045 - 2106 M 1230
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 116)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[255638, 255743, 255743, 255794],
                             [  1065,   1170,   1173,   1224]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 1065 1224 + gi|330443688|ref|NC_001145.3| 255638 255794 + 116 M 105 I 3 M 51
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_coding2genome_vulgar(self):
        """Test parsing exn_22_m_coding2genome_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_coding2genome_vulgar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments, check_operations=False)

    def check_vulgar(self, alignments, check_operations=True):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m coding2genome ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar no --showvulgar yes",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2151)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                # fmt: off
# flake8: noqa
                alignment.coordinates, numpy.array([[1318047, 1319274], [1228, 1]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"C"))
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: gi|296143771|ref|NM_001180731.1| 1228 1 - gi|330443520|ref|NC_001136.10| 1318047 1319274 + 2151 C 1227 1227
""",
            )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2106)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"C"))
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: gi|296143771|ref|NM_001180731.1| 0 1230 + gi|330443520|ref|NC_001136.10| 1319275 1318045 - 2106 C 1230 1230
""",
            )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 116)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[255638, 255743, 255743, 255794],
                             [  1065,   1170,   1173,   1224]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"CIC"))
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: gi|296143771|ref|NM_001180731.1| 1065 1224 + gi|330443688|ref|NC_001145.3| 255638 255794 + 116 C 105 105 G 3 0 C 51 51
""",
            )
        with self.assertRaises(StopIteration):
            next(alignments)


class Exonerate_dna2protein(unittest.TestCase):
    def test_exn_22_m_dna2protein_cigar(self):
        """Test parsing exn_22_m_dna2protein_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_dna2protein_cigar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)

    def check_cigar(self, alignments):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate --showcigar yes --showvulgar no --showalignment no nuc2.fa pro.fa",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "dna")
        self.assertEqual(alignment.target.id, "protein")
        self.assertEqual(alignment.score, 105)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[313, 344], [0, 93]]))
        )
        self.assertEqual(alignment.target.annotations["molecule_type"], "protein")
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: dna 0 93 + protein 313 344 . 105 M 93
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_dna2protein_vulgar(self):
        """Test parsing exn_22_m_dna2protein_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_dna2protein_vulgar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments, check_operations=False)

    def check_vulgar(self, alignments, check_operations=True):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate --showcigar no --showvulgar yes --showalignment no nuc2.fa pro.fa",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "dna")
        self.assertEqual(alignment.target.id, "protein")
        self.assertEqual(alignment.score, 105)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[313, 344], [0, 93]]))
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"M"))
        self.assertEqual(
            alignment.format("exonerate", "vulgar"),
            """\
vulgar: dna 0 93 + protein 313 344 . 105 M 93 31
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)


class Exonerate_genome2genome(unittest.TestCase):
    def test_exn_22_m_genome2genome_cigar(self):
        """Test parsing exn_22_o_vulgar_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_o_vulgar_cigar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 8)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments, check_operations=False)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 8)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)

    def check_cigar(self, alignments, check_operations=True):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m genome2genome ../intron.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showvulgar yes --showcigar yes",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sacCer3_dna")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2641)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[1319997, 1319971, 1319968, 1319468],
                             [    529,     503,     500,       0]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: sacCer3_dna 529 0 - gi|330443520|ref|NC_001136.10| 1319997 1319468 - 2641 M 26 M 3 M 500
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sacCer3_dna")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2641)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[1319997, 1319971, 1319968, 1319468],
                             [    529,     503,     500,       0]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"MCM"))
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: sacCer3_dna 529 0 - gi|330443520|ref|NC_001136.10| 1319997 1319468 - 2641 M 26 M 3 M 500
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sacCer3_dna")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2641)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[1319468, 1319558, 1319561, 1319997],
                             [      0,      90,      93,     529]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: sacCer3_dna 0 529 + gi|330443520|ref|NC_001136.10| 1319468 1319997 + 2641 M 90 M 3 M 436
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sacCer3_dna")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2641)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[1319468, 1319558, 1319561, 1319997],
                             [      0,      90,      93,     529]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"MCM"))
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: sacCer3_dna 0 529 + gi|330443520|ref|NC_001136.10| 1319468 1319997 + 2641 M 90 M 3 M 436
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sacCer3_dna")
        self.assertEqual(alignment.target.id, "gi|330443489|ref|NC_001135.5|")
        self.assertEqual(alignment.score, 267)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 23668,  23697,  32680,  32712,  32714,  32716,
                               32717,  32732,  42287,  42290,  42290,  42295,
                               42297,  42300,  42301,  42305,  42306,  42324,
                               42325,  97748,  97753,  97775,  97775,  97821,
                              115419, 115433, 115434, 115443, 115443, 115458,
                              115461, 115478, 115481, 115482, 115483, 115496,
                              115497, 115503, 115503, 115515, 115517, 115562,
                              115563, 115569],
                             [   491,    462,    462,    430,    430,    428,
                                 428,    413,    413,    410,    409,    404,
                                 404,    401,    401,    397,    397,    379,
                                 378,    378,    373,    351,    348,    302,
                                 302,    288,    288,    279,    278,    263,
                                 263,    246,    246,    245,    245,    232,
                                 232,    226,    225,    213,    213,    168,
                                 168,    162]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: sacCer3_dna 491 162 - gi|330443489|ref|NC_001135.5| 23668 115569 + 267 M 29 D 8983 M 32 D 2 M 2 D 1 M 15 D 9555 M 3 I 1 M 5 D 2 M 3 D 1 M 4 D 1 M 18 M 1 D 55423 M 5 M 22 I 3 M 46 D 17598 M 14 D 1 M 9 I 1 M 15 D 3 M 17 D 3 M 1 D 1 M 13 D 1 M 6 I 1 M 12 D 2 M 45 D 1 M 6
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sacCer3_dna")
        self.assertEqual(alignment.target.id, "gi|330443489|ref|NC_001135.5|")
        self.assertEqual(alignment.score, 267)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 23668,  23697,  23699,  32678,  32680,  32712,
                               32714,  32716,  32717,  32732,  32734,  42285,
                               42287,  42290,  42290,  42295,  42297,  42300,
                               42301,  42305,  42306,  42324,  42325,  42327,
                               97746,  97748,  97750,  97753,  97775,  97775,
                               97821,  97823, 115417, 115419, 115433, 115434,
                              115443, 115443, 115458, 115461, 115478, 115481,
                              115482, 115483, 115496, 115497, 115503, 115503,
                              115515, 115517, 115562, 115563, 115569],
                             [   491,    462,    462,    462,    462,    430,
                                 430,    428,    428,    413,    413,    413,
                                 413,    410,    409,    404,    404,    401,
                                 401,    397,    397,    379,    378,    378,
                                 378,    378,    376,    373,    351,    348,
                                 302,    302,    302,    302,    288,    288,
                                 279,    278,    263,    263,    246,    246,
                                 245,    245,    232,    232,    226,    225,
                                 213,    213,    168,    168,    162]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(
                alignment.operations,
                bytearray(b"M5N3MDMDM5N3MIMDMDMDMS5N3SCMIM5N3MDMIMDMDMDMDMIMDMDM"),
            )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: sacCer3_dna 491 162 - gi|330443489|ref|NC_001135.5| 23668 115569 + 267 M 29 D 2 D 8979 D 2 M 32 D 2 M 2 D 1 M 15 D 2 D 9551 D 2 M 3 I 1 M 5 D 2 M 3 D 1 M 4 D 1 M 18 D 1 I 1 D 2 D 55419 D 2 D 2 I 2 M 3 M 22 I 3 M 46 D 2 D 17594 D 2 M 14 D 1 M 9 I 1 M 15 D 3 M 17 D 3 M 1 D 1 M 13 D 1 M 6 I 1 M 12 D 2 M 45 D 1 M 6
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sacCer3_dna")
        self.assertEqual(alignment.target.id, "gi|330443667|ref|NC_001143.9|")
        self.assertEqual(alignment.score, 267)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[641760, 641729, 641729, 641725, 641725, 641706,
                              641703, 641694, 641693, 641687, 641687, 641680,
                              487436, 487436, 487389, 487362, 487362, 487358,
                              487358, 487355, 487355, 487351, 487351, 487342,
                              487341, 487325, 386209, 386209, 386184, 386183,
                              386168, 386168, 386159, 386159, 386157, 386157,
                              386143, 386123, 208677, 208664, 208662, 208661,
                              208637,  71940,  71940,  71934,  71934,  71933,
                               71933,  71932,  71932,  71931,  71931,  71930,
                               71930,  71929,  71929,  71928,  71928,  71927,
                               71927,  71926,  71926,  71925,  71925,  71924,
                               71924,  71923,  71923,  71921,  71921,  71919,
                               71919,  71905,  71905,  71883],
                             [   529,    498,    495,    491,    489,    470,
                                 470,    461,    461,    455,    454,    447,
                                 447,    390,    390,    363,    358,    354,
                                 353,    350,    347,    343,    342,    333,
                                 333,    317,    317,    286,    261,    261,
                                 246,    245,    236,    235,    233,    232,
                                 218,    198,    198,    185,    183,    183,
                                 159,    159,    152,    152,    151,    151,
                                 150,    150,    149,    149,    148,    148,
                                 146,    146,    145,    145,    144,    144,
                                 141,    141,    139,    139,    138,    138,
                                 137,    137,    135,    135,    133,    133,
                                 116,    102,    100,     78]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: sacCer3_dna 529 78 - gi|330443667|ref|NC_001143.9| 641760 71883 - 267 M 31 I 3 M 4 I 2 M 19 D 3 M 9 D 1 M 6 I 1 M 7 D 154244 I 57 D 47 M 27 I 5 M 4 I 1 M 3 I 3 M 4 I 1 M 9 D 1 M 16 D 101116 I 31 M 25 D 1 M 15 I 1 M 9 I 1 M 2 I 1 M 14 M 20 D 177446 M 13 M 2 D 1 M 24 D 136697 I 7 D 6 I 1 D 1 I 1 D 1 I 1 D 1 I 1 D 1 I 2 D 1 I 1 D 1 I 1 D 1 I 3 D 1 I 2 D 1 I 1 D 1 I 1 D 1 I 2 D 2 I 2 D 2 I 17 M 14 I 2 M 22
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sacCer3_dna")
        self.assertEqual(alignment.target.id, "gi|330443667|ref|NC_001143.9|")
        self.assertEqual(alignment.score, 267)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[641760, 641729, 641729, 641725, 641725, 641706,
                              641703, 641694, 641693, 641687, 641687, 641682,
                              641680, 487436, 487436, 487389, 487387, 487362,
                              487362, 487358, 487358, 487355, 487355, 487351,
                              487351, 487342, 487341, 487327, 487325, 386209,
                              386209, 386207, 386184, 386183, 386168, 386168,
                              386159, 386159, 386157, 386157, 386143, 386125,
                              386123, 386121, 208679, 208677, 208676, 208664,
                              208662, 208661, 208639, 208637,  71940,  71940,
                               71934,  71934,  71933,  71933,  71932,  71932,
                               71931,  71931,  71930,  71930,  71929,  71929,
                               71928,  71928,  71927,  71927,  71926,  71926,
                               71925,  71925,  71924,  71924,  71923,  71923,
                               71921,  71921,  71919,  71919,  71917,  71905,
                               71905,  71883],
                             [   529,    498,    495,    491,    489,    470,
                                 470,    461,    461,    455,    454,    449,
                                 447,    447,    390,    390,    388,    363,
                                 358,    354,    353,    350,    347,    343,
                                 342,    333,    333,    319,    317,    317,
                                 286,    284,    261,    261,    246,    245,
                                 236,    235,    233,    232,    218,    200,
                                 198,    198,    198,    198,    197,    185,
                                 183,    183,    161,    159,    159,    152,
                                 152,    151,    151,    150,    150,    149,
                                 149,    148,    148,    146,    146,    145,
                                 145,    144,    144,    141,    141,    139,
                                 139,    138,    138,    137,    137,    135,
                                 135,    133,    133,    116,    114,    102,
                                 100,     78]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(
                alignment.operations,
                bytearray(
                    b"MIMIMDMDMIM5NNN3MIMIMIMIMDM5NN3MDMIMIMIMCS5N3SCMDM5NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN3MIM"
                ),
            )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: sacCer3_dna 529 78 - gi|330443667|ref|NC_001143.9| 641760 71883 - 267 M 31 I 3 M 4 I 2 M 19 D 3 M 9 D 1 M 6 I 1 M 5 M 2 D 154244 I 57 D 47 M 2 M 25 I 5 M 4 I 1 M 3 I 3 M 4 I 1 M 9 D 1 M 14 M 2 D 101116 I 31 M 2 M 23 D 1 M 15 I 1 M 9 I 1 M 2 I 1 M 14 M 18 D 2 I 2 D 2 D 177442 D 2 D 1 I 1 M 12 M 2 D 1 M 22 M 2 D 136697 I 7 D 6 I 1 D 1 I 1 D 1 I 1 D 1 I 1 D 1 I 2 D 1 I 1 D 1 I 1 D 1 I 3 D 1 I 2 D 1 I 1 D 1 I 1 D 1 I 2 D 2 I 2 D 2 I 17 M 2 M 12 I 2 M 22
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_genome2genome_vulgar(self):
        """Test parsing exn_22_o_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_o_vulgar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 4)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments)

    def check_vulgar(self, alignments, check_operations=True):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m genome2genome ../intron.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showvulgar yes --showcigar no",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sacCer3_dna")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2641)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[1319997, 1319971, 1319968, 1319468],
                             [    529,     503,     500,       0]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"MCM"))
        self.assertEqual(
            alignment.format("exonerate", "vulgar"),
            """\
vulgar: sacCer3_dna 529 0 - gi|330443520|ref|NC_001136.10| 1319997 1319468 - 2641 M 26 26 C 3 3 M 500 500
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sacCer3_dna")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2641)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[1319468, 1319558, 1319561, 1319997], [0, 90, 93, 529]]),
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"MCM"))
        self.assertEqual(
            alignment.format("exonerate", "vulgar"),
            """\
vulgar: sacCer3_dna 0 529 + gi|330443520|ref|NC_001136.10| 1319468 1319997 + 2641 M 90 90 C 3 3 M 436 436
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sacCer3_dna")
        self.assertEqual(alignment.target.id, "gi|330443489|ref|NC_001135.5|")
        self.assertEqual(alignment.score, 267)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 23668,  23697,  23699,  32678,  32680,  32712,
                               32714,  32716,  32717,  32732,  32734,  42285,
                               42287,  42290,  42290,  42295,  42297,  42300,
                               42301,  42305,  42306,  42324,  42325,  42327,
                               97746,  97748,  97750,  97753,  97775,  97775,
                               97821,  97823, 115417, 115419, 115433, 115434,
                              115443, 115443, 115458, 115461, 115478, 115481,
                              115482, 115483, 115496, 115497, 115503, 115503,
                              115515, 115517, 115562, 115563, 115569],
                             [   491,    462,    462,    462,    462,    430,
                                 430,    428,    428,    413,    413,    413,
                                 413,    410,    409,    404,    404,    401,
                                 401,    397,    397,    379,    378,    378,
                                 378,    378,    376,    373,    351,    348,
                                 302,    302,    302,    302,    288,    288,
                                 279,    278,    263,    263,    246,    246,
                                 245,    245,    232,    232,    226,    225,
                                 213,    213,    168,    168,    162]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(
                alignment.operations,
                bytearray(b"M5N3MDMDM5N3MIMDMDMDMS5N3SCMIM5N3MDMIMDMDMDMDMIMDMDM"),
            )
        self.assertEqual(
            alignment.format("exonerate", "vulgar"),
            """\
vulgar: sacCer3_dna 491 162 - gi|330443489|ref|NC_001135.5| 23668 115569 + 267 M 29 29 5 0 2 I 0 8979 3 0 2 M 32 32 G 0 2 M 2 2 G 0 1 M 15 15 5 0 2 I 0 9551 3 0 2 M 3 3 G 1 0 M 5 5 G 0 2 M 3 3 G 0 1 M 4 4 G 0 1 M 18 18 S 1 1 5 0 2 I 0 55419 3 0 2 S 2 2 C 3 3 M 22 22 G 3 0 M 46 46 5 0 2 I 0 17594 3 0 2 M 14 14 G 0 1 M 9 9 G 1 0 M 15 15 G 0 3 M 17 17 G 0 3 M 1 1 G 0 1 M 13 13 G 0 1 M 6 6 G 1 0 M 12 12 G 0 2 M 45 45 G 0 1 M 6 6
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sacCer3_dna")
        self.assertEqual(alignment.target.id, "gi|330443667|ref|NC_001143.9|")
        self.assertEqual(alignment.score, 267)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[641760, 641729, 641729, 641725, 641725, 641706,
                              641703, 641694, 641693, 641687, 641687, 641682,
                              641680, 487436, 487436, 487389, 487387, 487362,
                              487362, 487358, 487358, 487355, 487355, 487351,
                              487351, 487342, 487341, 487327, 487325, 386209,
                              386209, 386207, 386184, 386183, 386168, 386168,
                              386159, 386159, 386157, 386157, 386143, 386125,
                              386123, 386121, 208679, 208677, 208676, 208664,
                              208662, 208661, 208639, 208637,  71940,  71940,
                               71934,  71934,  71933,  71933,  71932,  71932,
                               71931,  71931,  71930,  71930,  71929,  71929,
                               71928,  71928,  71927,  71927,  71926,  71926,
                               71925,  71925,  71924,  71924,  71923,  71923,
                               71921,  71921,  71919,  71919,  71917,  71905,
                               71905,  71883],
                             [   529,    498,    495,    491,    489,    470,
                                 470,    461,    461,    455,    454,    449,
                                 447,    447,    390,    390,    388,    363,
                                 358,    354,    353,    350,    347,    343,
                                 342,    333,    333,    319,    317,    317,
                                 286,    284,    261,    261,    246,    245,
                                 236,    235,    233,    232,    218,    200,
                                 198,    198,    198,    198,    197,    185,
                                 183,    183,    161,    159,    159,    152,
                                 152,    151,    151,    150,    150,    149,
                                 149,    148,    148,    146,    146,    145,
                                 145,    144,    144,    141,    141,    139,
                                 139,    138,    138,    137,    137,    135,
                                 135,    133,    133,    116,    114,    102,
                                 100,     78]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(
                alignment.operations,
                bytearray(
                    b"MIMIMDMDMIM5NNN3MIMIMIMIMDM5NN3MDMIMIMIMCS5N3SCMDM5NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN3MIM"
                ),
            )
        self.assertEqual(
            alignment.format("exonerate", "vulgar"),
            """\
vulgar: sacCer3_dna 529 78 - gi|330443667|ref|NC_001143.9| 641760 71883 - 267 M 31 31 G 3 0 M 4 4 G 2 0 M 19 19 G 0 3 M 9 9 G 0 1 M 6 6 G 1 0 M 5 5 5 2 2 I 0 154244 I 57 0 I 0 47 3 2 2 M 25 25 G 5 0 M 4 4 G 1 0 M 3 3 G 3 0 M 4 4 G 1 0 M 9 9 G 0 1 M 14 14 5 2 2 I 0 101116 I 31 0 3 2 2 M 23 23 G 0 1 M 15 15 G 1 0 M 9 9 G 1 0 M 2 2 G 1 0 M 14 14 C 18 18 S 2 2 5 0 2 I 0 177442 3 0 2 S 1 1 C 12 12 M 2 2 G 0 1 M 22 22 5 2 2 I 0 136697 I 7 0 I 0 6 I 1 0 I 0 1 I 1 0 I 0 1 I 1 0 I 0 1 I 1 0 I 0 1 I 2 0 I 0 1 I 1 0 I 0 1 I 1 0 I 0 1 I 3 0 I 0 1 I 2 0 I 0 1 I 1 0 I 0 1 I 1 0 I 0 1 I 2 0 I 0 2 I 2 0 I 0 2 I 17 0 3 2 2 M 12 12 G 2 0 M 22 22
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)


class Exonerate_ungapped(unittest.TestCase):
    def test_exn_22_m_ungapped_cigar(self):
        """Test parsing exn_22_m_ungapped_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_ungapped_cigar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)

    def check_cigar(self, alignments):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m ungapped ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar yes --showvulgar no",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 6150)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 0 1230 + gi|330443520|ref|NC_001136.10| 1319275 1318045 - 6150 M 1230
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 233)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[254031, 254146], [121, 236]])
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 121 236 + gi|330443688|ref|NC_001145.3| 254031 254146 + 233 M 115
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 151)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[255671, 255739], [1098, 1166]])
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 1098 1166 + gi|330443688|ref|NC_001145.3| 255671 255739 + 151 M 68
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_ungapped_vulgar(self):
        """Test parsing exn_22_m_ungapped_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_ungapped_vulgar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments, check_operations=False)

    def check_vulgar(self, alignments, check_operations=True):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m ungapped ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar no --showvulgar yes",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 6150)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"M"))
        self.assertEqual(
            alignment.format("exonerate", "vulgar"),
            """\
vulgar: gi|296143771|ref|NM_001180731.1| 0 1230 + gi|330443520|ref|NC_001136.10| 1319275 1318045 - 6150 M 1230 1230
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 233)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[254031, 254146], [121, 236]])
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"M"))
        self.assertEqual(
            alignment.format("exonerate", "vulgar"),
            """\
vulgar: gi|296143771|ref|NM_001180731.1| 121 236 + gi|330443688|ref|NC_001145.3| 254031 254146 + 233 M 115 115
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 151)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[255671, 255739], [1098, 1166]])
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"M"))
        self.assertEqual(
            alignment.format("exonerate", "vulgar"),
            """\
vulgar: gi|296143771|ref|NM_001180731.1| 1098 1166 + gi|330443688|ref|NC_001145.3| 255671 255739 + 151 M 68 68
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)


class Exonerate_ungapped_trans(unittest.TestCase):
    def test_exn_22_m_ungapped_trans_cigar(self):
        """Test parsing exn_22_m_ungapped_trans_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_ungapped_trans_cigar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)

    def check_cigar(self, alignments):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m ungapped:trans ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar yes --showvulgar no",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2151)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1318047, 1319274], [1228, 1]])
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 1228 1 - gi|330443520|ref|NC_001136.10| 1318047 1319274 + 2151 M 1227
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2106)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 0 1230 + gi|330443520|ref|NC_001136.10| 1319275 1318045 - 2106 M 1230
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2072)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1318045, 1319275], [1230, 0]])
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 1230 0 - gi|330443520|ref|NC_001136.10| 1318045 1319275 + 2072 M 1230
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_ungapped_trans_vulgar(self):
        """Test parsing exn_22_m_ungapped_trans_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_ungapped_trans_vulgar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments, check_operations=False)

    def check_vulgar(self, alignments, check_operations=True):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m ungapped:trans ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar no --showvulgar yes",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2151)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1318047, 1319274], [1228, 1]])
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"C"))
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: gi|296143771|ref|NM_001180731.1| 1228 1 - gi|330443520|ref|NC_001136.10| 1318047 1319274 + 2151 C 1227 1227
""",
            )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2106)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"C"))
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: gi|296143771|ref|NM_001180731.1| 0 1230 + gi|330443520|ref|NC_001136.10| 1319275 1318045 - 2106 C 1230 1230
""",
            )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2072)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1318045, 1319275], [1230, 0]])
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"C"))
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: gi|296143771|ref|NM_001180731.1| 1230 0 - gi|330443520|ref|NC_001136.10| 1318045 1319275 + 2072 C 1230 1230
""",
            )
        with self.assertRaises(StopIteration):
            next(alignments)


class Exonerate_ner(unittest.TestCase):
    def test_exn_22_m_ner_cigar(self):
        """Test parsing exonerate output (exn_22_m_ner_cigar.exn)."""
        exn_file = os.path.join("Exonerate", "exn_22_m_ner_cigar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)

    def check_cigar(self, alignments):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m ner ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar yes --showvulgar no",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 6150)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 0 1230 + gi|330443520|ref|NC_001136.10| 1319275 1318045 - 6150 M 1230
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443681|ref|NC_001144.5|")
        self.assertEqual(alignment.score, 502)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[297910, 297918, 297946, 297946, 297958, 297970,
                              297970, 297984, 297992, 297992, 297997, 297997,
                              298005, 298016, 298016, 298019, 298023, 298024,
                              298031, 298032, 298039, 298049, 298049, 298056,
                              298065, 298065, 298084, 298095, 298095, 298108,
                              298117, 298123, 298123, 298133, 298135, 298139,
                              298140, 298150, 298150, 318437, 318446, 318452,
                              318452, 318454, 318461, 318476, 318476, 318488,
                              318488, 318498, 318503, 318503, 318523, 318523,
                              318529, 318532, 318532, 318551, 318551, 318558,
                              318558, 318559, 318559, 318560, 318560, 318562,
                              318562, 318563, 318563, 318564, 318564, 318565,
                              318565, 318571, 318611, 318611, 318625, 318630,
                              318630, 318638, 318639, 318639, 318646, 318646,
                              318655, 318659, 318659, 318666, 318681, 318681,
                              318688, 318691, 318700, 318708, 318708, 318713,
                              318713, 318722, 318725, 318725, 318735, 318742,
                              318749, 318755, 318755, 318766, 318767, 318773,
                              318787, 318787, 318793, 318793, 318808, 318808,
                              318822, 318834, 318834, 318840, 318851, 318851,
                              318871, 318889, 318889, 318893, 318893, 318898,
                              318913, 318913, 318920, 318922, 318922, 318924,
                              318934, 318944, 318944, 318955, 318959, 318959,
                              318973, 318974, 318974, 318994],
                             [   110,    118,    118,    148,    160,    160,
                                 169,    183,    183,    184,    189,    190,
                                 198,    198,    227,    227,    231,    231,
                                 238,    238,    245,    245,    255,    262,
                                 262,    266,    285,    285,    296,    309,
                                 309,    315,    316,    326,    326,    330,
                                 330,    340,    708,    708,    717,    717,
                                 728,    728,    735,    735,    744,    756,
                                 758,    768,    768,    778,    798,    799,
                                 805,    805,    806,    806,    807,    807,
                                 808,    808,    809,    809,    812,    812,
                                 813,    813,    814,    814,    815,    815,
                                 816,    822,    822,    832,    846,    846,
                                 889,    897,    897,    915,    922,    923,
                                 932,    932,    940,    947,    947,    953,
                                 960,    960,    969,    969,    979,    984,
                                 985,    994,    994,   1006,   1016,   1016,
                                1023,   1023,   1024,   1035,   1035,   1041,
                                1041,   1047,   1053,   1054,   1069,   1070,
                                1084,   1084,   1092,   1098,   1098,   1099,
                                1119,   1119,   1121,   1125,   1126,   1131,
                                1131,   1133,   1140,   1140,   1156,   1156,
                                1166,   1166,   1167,   1178,   1178,   1187,
                                1201,   1201,   1210,   1230]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 110 1230 + gi|330443681|ref|NC_001144.5| 297910 318994 + 502 M 8 D 28 I 30 M 12 D 12 I 9 M 14 D 8 I 1 M 5 I 1 M 8 D 11 I 29 D 3 M 4 D 1 M 7 D 1 M 7 D 10 I 10 M 7 D 9 I 4 M 19 D 11 I 11 M 13 D 9 M 6 I 1 M 10 D 2 M 4 D 1 M 10 I 368 D 20287 M 9 D 6 I 11 D 2 M 7 D 15 I 9 M 12 I 2 M 10 D 5 I 10 M 20 I 1 M 6 D 3 I 1 D 19 I 1 D 7 I 1 D 1 I 1 D 1 I 3 D 2 I 1 D 1 I 1 D 1 I 1 D 1 I 1 M 6 D 40 I 10 M 14 D 5 I 43 M 8 D 1 I 18 M 7 I 1 M 9 D 4 I 8 M 7 D 15 I 6 M 7 D 3 M 9 D 8 I 10 M 5 I 1 M 9 D 3 I 12 M 10 D 7 M 7 D 6 I 1 M 11 D 1 M 6 D 14 I 6 M 6 I 1 M 15 I 1 M 14 D 12 I 8 M 6 D 11 I 1 M 20 D 18 I 2 M 4 I 1 M 5 D 15 I 2 M 7 D 2 I 16 D 2 M 10 D 10 I 1 M 11 D 4 I 9 M 14 D 1 I 9 M 20
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 440)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[183946, 183952, 183960, 183960, 183977, 183988,
                              183995, 183995, 184002, 184031, 184039, 184039,
                              184044, 184055, 184059, 184059, 184060, 184060,
                              184063, 184063, 184064, 184064, 184066, 184066,
                              184081, 184083, 184083, 184092, 184108, 184111,
                              184118, 184119, 184125, 184126, 184132, 184132,
                              184148, 184169, 184169, 184179, 184190, 184190,
                              184194, 184196, 184200, 184200, 184210, 184211,
                              184211, 184214, 184214, 184225, 184239, 184239,
                              184252, 184271, 184271, 184280, 184280, 184289,
                              184299, 184299, 184305, 184311, 184311, 184314,
                              184314, 184315, 184315, 184337, 184339, 184339,
                              184346, 184365, 184366, 184366, 184378, 184390,
                              184390, 184393, 184399, 184399, 184410, 184411,
                              184415, 184417, 184428, 184428, 184433, 184439,
                              184439, 184450, 184463, 184463, 184473, 184480,
                              184480, 184481, 184481, 184483, 184483, 184484,
                              184484, 184485, 184485, 184495, 184495, 184514,
                              184515, 184515, 184525, 184539, 184539, 184541,
                              184541, 184542, 184542, 184552, 184554, 184563,
                              184565, 184565, 184573, 184576, 184576, 184578,
                              184603],
                             [   509,    515,    515,    537,    537,    548,
                                 548,    567,    567,    596,    596,    607,
                                 607,    618,    618,    619,    619,    621,
                                 621,    623,    623,    626,    626,    636,
                                 651,    651,    667,    667,    683,    683,
                                 690,    690,    696,    696,    702,    703,
                                 719,    719,    737,    747,    747,    748,
                                 752,    752,    756,    757,    767,    767,
                                 777,    780,    781,    792,    792,    797,
                                 810,    810,    815,    824,    825,    834,
                                 834,    843,    849,    849,    852,    852,
                                 853,    853,    861,    883,    883,    902,
                                 902,    921,    921,    927,    939,    939,
                                 957,    957,    963,    964,    975,    975,
                                 979,    979,    990,    991,    996,    996,
                                1010,   1021,   1021,   1022,   1032,   1032,
                                1034,   1034,   1035,   1035,   1037,   1037,
                                1039,   1039,   1050,   1060,   1061,   1080,
                                1080,   1099,   1109,   1109,   1110,   1110,
                                1112,   1112,   1122,   1132,   1132,   1141,
                                1141,   1147,   1155,   1155,   1167,   1167,
                                1192]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 509 1192 + gi|330443520|ref|NC_001136.10| 183946 184603 + 440 M 6 D 8 I 22 D 17 M 11 D 7 I 19 D 7 M 29 D 8 I 11 D 5 M 11 D 4 I 1 D 1 I 2 D 3 I 2 D 1 I 3 D 2 I 10 M 15 D 2 I 16 D 9 M 16 D 3 M 7 D 1 M 6 D 1 M 6 I 1 M 16 D 21 I 18 M 10 D 11 I 1 M 4 D 2 M 4 I 1 M 10 D 1 I 10 M 3 I 1 M 11 D 14 I 5 M 13 D 19 I 5 M 9 I 1 M 9 D 10 I 9 M 6 D 6 I 3 D 3 I 1 D 1 I 8 M 22 D 2 I 19 D 7 M 19 D 1 I 6 M 12 D 12 I 18 D 3 M 6 I 1 M 11 D 1 M 4 D 2 M 11 I 1 M 5 D 6 I 14 M 11 D 13 I 1 M 10 D 7 I 2 D 1 I 1 D 2 I 2 D 1 I 2 D 1 I 11 M 10 I 1 M 19 D 1 I 19 M 10 D 14 I 1 D 2 I 2 D 1 I 10 M 10 D 2 M 9 D 2 I 6 M 8 D 3 I 12 D 2 M 25
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_ner_vulgar(self):
        """Test parsing exonerate output (exn_22_m_ner_vulgar.exn)."""
        exn_file = os.path.join("Exonerate", "exn_22_m_ner_vulgar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments, check_operations=False)

    def check_vulgar(self, alignments, check_operations=True):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m ner ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar no --showvulgar yes",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 6150)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"M"))
        self.assertEqual(
            alignment.format("exonerate", "vulgar"),
            """\
vulgar: gi|296143771|ref|NM_001180731.1| 0 1230 + gi|330443520|ref|NC_001136.10| 1319275 1318045 - 6150 M 1230 1230
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443681|ref|NC_001144.5|")
        self.assertEqual(alignment.score, 502)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[297910, 297917, 297946, 297946, 297957, 297970,
                              297970, 297983, 297992, 297992, 297997, 297997,
                              298004, 298019, 298019, 298023, 298024, 298031,
                              298032, 298038, 298049, 298049, 298055, 298065,
                              298065, 298083, 298095, 298095, 298107, 298117,
                              298117, 298123, 298123, 298133, 298135, 298139,
                              298140, 298149, 318437, 318437, 318445, 318454,
                              318454, 318460, 318476, 318476, 318488, 318488,
                              318497, 318503, 318503, 318523, 318523, 318528,
                              318565, 318565, 318570, 318611, 318611, 318624,
                              318630, 318630, 318637, 318639, 318639, 318646,
                              318646, 318654, 318659, 318659, 318665, 318681,
                              318681, 318687, 318691, 318691, 318699, 318708,
                              318708, 318713, 318713, 318721, 318725, 318725,
                              318734, 318742, 318742, 318748, 318755, 318755,
                              318766, 318767, 318772, 318787, 318787, 318793,
                              318793, 318808, 318808, 318821, 318834, 318834,
                              318839, 318851, 318851, 318870, 318889, 318889,
                              318893, 318893, 318897, 318913, 318913, 318919,
                              318924, 318924, 318933, 318944, 318944, 318954,
                              318959, 318959, 318972, 318974, 318974, 318994],
                             [   110,    117,    117,    148,    159,    159,
                                 169,    182,    182,    184,    189,    190,
                                 197,    197,    227,    231,    231,    238,
                                 238,    244,    244,    255,    261,    261,
                                 266,    284,    284,    296,    308,    308,
                                 309,    315,    316,    326,    326,    330,
                                 330,    339,    339,    708,    716,    716,
                                 728,    734,    734,    744,    756,    758,
                                 767,    767,    778,    798,    799,    804,
                                 804,    816,    821,    821,    832,    845,
                                 845,    889,    896,    896,    915,    922,
                                 923,    931,    931,    940,    946,    946,
                                 953,    959,    959,    960,    968,    968,
                                 979,    984,    985,    993,    993,   1006,
                                1015,   1015,   1016,   1022,   1022,   1024,
                                1035,   1035,   1040,   1040,   1047,   1053,
                                1054,   1069,   1070,   1083,   1083,   1092,
                                1097,   1097,   1099,   1118,   1118,   1121,
                                1125,   1126,   1130,   1130,   1133,   1139,
                                1139,   1156,   1165,   1165,   1167,   1177,
                                1177,   1187,   1200,   1200,   1210,   1230]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(
                alignment.operations,
                bytearray(
                    b"MUUMUUMUUMIMUUMDMDMUUMUUMUUMUUMIMDMDMUUMUUMUUMIMUUMIMUUMUUMUUMUUMIMUUMUUMUUMUUMIMUUMUUMUUMDMUUMIMIMUUMUUMUUMIMUUMUUMUUMUUMUUM",
                ),
            )
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: gi|296143771|ref|NM_001180731.1| 110 1230 + gi|330443681|ref|NC_001144.5| 297910 318994 + 502 M 7 7 N 31 29 M 11 11 N 10 13 M 13 13 N 2 9 M 5 5 G 1 0 M 7 7 N 30 15 M 4 4 G 0 1 M 7 7 G 0 1 M 6 6 N 11 11 M 6 6 N 5 10 M 18 18 N 12 12 M 12 12 N 1 10 M 6 6 G 1 0 M 10 10 G 0 2 M 4 4 G 0 1 M 9 9 N 369 20288 M 8 8 N 12 9 M 6 6 N 10 16 M 12 12 G 2 0 M 9 9 N 11 6 M 20 20 G 1 0 M 5 5 N 12 37 M 5 5 N 11 41 M 13 13 N 44 6 M 7 7 N 19 2 M 7 7 G 1 0 M 8 8 N 9 5 M 6 6 N 7 16 M 6 6 N 1 4 M 8 8 N 11 9 M 5 5 G 1 0 M 8 8 N 13 4 M 9 9 N 1 8 M 6 6 N 2 7 M 11 11 G 0 1 M 5 5 N 7 15 M 6 6 G 1 0 M 15 15 G 1 0 M 13 13 N 9 13 M 5 5 N 2 12 M 19 19 N 3 19 M 4 4 G 1 0 M 4 4 N 3 16 M 6 6 N 17 5 M 9 9 N 2 11 M 10 10 N 10 5 M 13 13 N 10 2 M 20 20
""",
            )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 440)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[183946, 183951, 183977, 183977, 183987, 184002,
                              184002, 184030, 184044, 184044, 184054, 184066,
                              184066, 184080, 184092, 184092, 184107, 184111,
                              184111, 184118, 184119, 184125, 184126, 184132,
                              184132, 184147, 184169, 184169, 184178, 184190,
                              184190, 184194, 184196, 184200, 184200, 184209,
                              184211, 184211, 184214, 184214, 184224, 184239,
                              184239, 184251, 184271, 184271, 184280, 184280,
                              184288, 184299, 184299, 184304, 184315, 184315,
                              184336, 184346, 184346, 184364, 184366, 184366,
                              184377, 184393, 184393, 184399, 184399, 184410,
                              184411, 184415, 184417, 184428, 184428, 184432,
                              184439, 184439, 184449, 184463, 184463, 184472,
                              184485, 184485, 184495, 184495, 184513, 184515,
                              184515, 184524, 184542, 184542, 184552, 184554,
                              184562, 184565, 184565, 184572, 184578, 184578,
                              184603],
                             [   509,    514,    514,    537,    547,    547,
                                 567,    595,    595,    607,    617,    617,
                                 636,    650,    650,    667,    682,    682,
                                 683,    690,    690,    696,    696,    702,
                                 703,    718,    718,    737,    746,    746,
                                 748,    752,    752,    756,    757,    766,
                                 766,    777,    780,    781,    791,    791,
                                 797,    809,    809,    815,    824,    825,
                                 833,    833,    843,    848,    848,    861,
                                 882,    882,    902,    920,    920,    927,
                                 938,    938,    957,    963,    964,    975,
                                 975,    979,    979,    990,    991,    995,
                                 995,   1010,   1020,   1020,   1022,   1031,
                                1031,   1050,   1060,   1061,   1079,   1079,
                                1099,   1108,   1108,   1122,   1132,   1132,
                                1140,   1140,   1147,   1154,   1154,   1167,
                                1192]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(
                alignment.operations,
                bytearray(
                    b"MUUMUUMUUMUUMUUMUUMDMDMIMUUMUUMDMIMUUMIMUUMUUMIMUUMUUMUUMUUMUUMIMDMDMIMUUMUUMUUMIMUUMUUMDMUUMUUM"
                ),
            )
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: gi|296143771|ref|NM_001180731.1| 509 1192 + gi|330443520|ref|NC_001136.10| 183946 184603 + 440 M 5 5 N 23 26 M 10 10 N 20 15 M 28 28 N 12 14 M 10 10 N 19 12 M 14 14 N 17 12 M 15 15 N 1 4 M 7 7 G 0 1 M 6 6 G 0 1 M 6 6 G 1 0 M 15 15 N 19 22 M 9 9 N 2 12 M 4 4 G 0 2 M 4 4 G 1 0 M 9 9 N 11 2 M 3 3 G 1 0 M 10 10 N 6 15 M 12 12 N 6 20 M 9 9 G 1 0 M 8 8 N 10 11 M 5 5 N 13 11 M 21 21 N 20 10 M 18 18 N 7 2 M 11 11 N 19 16 M 6 6 G 1 0 M 11 11 G 0 1 M 4 4 G 0 2 M 11 11 G 1 0 M 4 4 N 15 7 M 10 10 N 2 14 M 9 9 N 19 13 M 10 10 G 1 0 M 18 18 N 20 2 M 9 9 N 14 18 M 10 10 G 0 2 M 8 8 N 7 3 M 7 7 N 13 6 M 25 25
""",
            )
        with self.assertRaises(StopIteration):
            next(alignments)


class Exonerate_multiple(unittest.TestCase):
    def test_exn_22_q_multiple_cigar(self):
        """Test parsing exn_22_q_multiple_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_q_multiple_cigar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 6)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 6)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)

    def check_cigar(self, alignments):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m est2genome comb.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar yes --showvulgar no",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.target.id, "gi|330443482|ref|NC_001134.8|")
        self.assertEqual(alignment.score, 4485)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[560077, 560974], [0, 897]])
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296142823|ref|NM_001178508.1| 0 897 + gi|330443482|ref|NC_001134.8| 560077 560974 + 4485 M 897
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.target.id, "gi|330443753|ref|NC_001148.4|")
        self.assertEqual(alignment.score, 941)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[492933, 492929, 492929, 492908, 492908, 492857,
                              492857, 492851, 492851, 492840, 492840, 492839,
                              492839, 492813, 492813, 492711, 492711, 492710,
                              492710, 492709, 492709, 492654, 492651, 492646,
                              492646, 492639, 492638, 492634, 492633, 492624,
                              492624, 492621, 492620, 492618, 492616, 492614,
                              492611, 492608, 492607, 492585, 492584, 492575,
                              492575, 492507, 492506, 492503, 492503, 492329,
                              492329, 492321, 492320, 492166, 492163, 492155,
                              492154, 492150, 492148, 492139, 492137, 492136,
                              492135, 492130, 492124, 492071, 492070, 492056,
                              492056, 492033],
                             [     2,      6,      8,     29,     30,     81,
                                  83,     89,     92,    103,    106,    107,
                                 108,    134,    137,    239,    240,    241,
                                 242,    243,    244,    299,    299,    304,
                                 306,    313,    313,    317,    317,    326,
                                 327,    330,    330,    332,    332,    334,
                                 334,    337,    337,    359,    359,    368,
                                 369,    437,    437,    440,    441,    615,
                                 616,    624,    624,    778,    778,    786,
                                 786,    790,    790,    799,    799,    800,
                                 800,    805,    805,    858,    858,    872,
                                 873,    896]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296142823|ref|NM_001178508.1| 2 896 + gi|330443753|ref|NC_001148.4| 492933 492033 - 941 M 4 I 2 M 21 I 1 M 51 I 2 M 6 I 3 M 11 I 3 M 1 I 1 M 26 I 3 M 102 I 1 M 1 I 1 M 1 I 1 M 55 D 3 M 5 I 2 M 7 D 1 M 4 D 1 M 9 I 1 M 3 D 1 M 2 D 2 M 2 D 3 M 3 D 1 M 22 D 1 M 9 I 1 M 68 D 1 M 3 I 1 M 174 I 1 M 8 D 1 M 154 D 3 M 8 D 1 M 4 D 2 M 9 D 2 M 1 D 1 M 5 D 6 M 53 D 1 M 14 I 1 M 23
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 651)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[267809, 267843, 267843, 267851, 267851, 267852,
                              267852, 267859, 267863, 267867, 267867, 267870,
                              267870, 267880, 267880, 267908, 267908, 267942,
                              267942, 267948, 267950, 267960, 267961, 267975,
                              267976, 267991, 267991, 268003, 268003, 268004,
                              268004, 268007, 268007, 268035, 268035, 268038,
                              268038, 268059, 268059, 268066, 268070, 268072,
                              268074, 268080, 268083, 268096, 268098, 268102,
                              268103, 268114, 268114, 268120, 268120, 268132,
                              268133, 268185, 268185, 268192, 268193, 268248,
                              268249, 268256, 268256, 268269, 268271, 268371,
                              268372, 268378, 268378, 268418, 268418, 268424,
                              268426, 268448, 300686, 300698, 300699, 300717],
                             [    34,     68,     69,     77,     78,     79,
                                  81,     88,     88,     92,     94,     97,
                                 100,    110,    113,    141,    142,    176,
                                 177,    183,    183,    193,    193,    207,
                                 207,    222,    223,    235,    236,    237,
                                 238,    241,    244,    272,    280,    283,
                                 285,    306,    308,    315,    315,    317,
                                 317,    323,    323,    336,    336,    340,
                                 340,    351,    356,    362,    364,    376,
                                 376,    428,    429,    436,    436,    491,
                                 491,    498,    501,    514,    514,    614,
                                 614,    620,    621,    661,    663,    669,
                                 669,    691,    691,    703,    703,    721]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296142823|ref|NM_001178508.1| 34 721 + gi|330443520|ref|NC_001136.10| 267809 300717 + 651 M 34 I 1 M 8 I 1 M 1 I 2 M 7 D 4 M 4 I 2 M 3 I 3 M 10 I 3 M 28 I 1 M 34 I 1 M 6 D 2 M 10 D 1 M 14 D 1 M 15 I 1 M 12 I 1 M 1 I 1 M 3 I 3 M 28 I 8 M 3 I 2 M 21 I 2 M 7 D 4 M 2 D 2 M 6 D 3 M 13 D 2 M 4 D 1 M 11 I 5 M 6 I 2 M 12 D 1 M 52 I 1 M 7 D 1 M 55 D 1 M 7 I 3 M 13 D 2 M 100 D 1 M 6 I 1 M 40 I 2 M 6 D 2 M 22 D 32238 M 12 D 1 M 18
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 6150)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 0 1230 + gi|330443520|ref|NC_001136.10| 1319275 1318045 - 6150 M 1230
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 439)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 85010,  85021,  85021,  85036,  85036,  85040,
                               85040,  85041,  85041,  85049,  85049,  85066,
                              253974, 253978, 253979, 253987, 253987, 253990,
                              253990, 254023, 254024, 254031, 254033, 254135,
                              350959, 350973, 350975, 350985, 350985, 350990,
                              350992, 351002, 351002, 351006, 351007, 351027,
                              351027, 351042, 351043, 351048, 351048, 351052,
                              473170, 473190, 473195, 473201],
                             [     0,     11,     12,     27,     29,     33,
                                  34,     35,     36,     44,     48,     65,
                                  65,     69,     69,     77,     78,     81,
                                  83,    116,    116,    123,    123,    225,
                                 225,    239,    239,    249,    251,    256,
                                 256,    266,    268,    272,    272,    292,
                                 293,    308,    308,    313,    316,    320,
                                 320,    340,    340,    346]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 0 346 + gi|330443688|ref|NC_001145.3| 85010 473201 + 439 M 11 I 1 M 15 I 2 M 4 I 1 M 1 I 1 M 8 I 4 M 17 D 168908 M 4 D 1 M 8 I 1 M 3 I 2 M 33 D 1 M 7 D 2 M 102 D 96824 M 14 D 2 M 10 I 2 M 5 D 2 M 10 I 2 M 4 D 1 M 20 I 1 M 15 D 1 M 5 I 3 M 4 D 122118 M 20 D 5 M 6
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 263)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[130198, 130184, 130183, 130179, 130179, 130154,
                              130153, 130144, 130138, 130096, 130096, 130080,
                              130078, 130071, 130070, 130067, 130067, 130044,
                              130044, 130038, 120681, 120680, 120680, 120669,
                              120668, 120656, 120656, 120647, 120646, 120636,
                              120636, 120618, 120617, 120612,  11487,  11471,
                               11471,  11467,  11467,  11456,  11456,  11448,
                               11448,  11426,  11424,  11420,  11418,  11384,
                               11383,  11380,  11379,  11338],
                             [    25,     39,     39,     43,     45,     70,
                                  70,     79,     79,    121,    123,    139,
                                 139,    146,    146,    149,    151,    174,
                                 177,    183,    183,    184,    185,    196,
                                 196,    208,    209,    218,    218,    228,
                                 229,    247,    247,    252,    252,    268,
                                 272,    276,    277,    288,    293,    301,
                                 302,    324,    324,    328,    328,    362,
                                 362,    365,    365,    406]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 25 406 + gi|330443688|ref|NC_001145.3| 130198 11338 - 263 M 14 D 1 M 4 I 2 M 25 D 1 M 9 D 6 M 42 I 2 M 16 D 2 M 7 D 1 M 3 I 2 M 23 I 3 M 6 D 9357 M 1 I 1 M 11 D 1 M 12 I 1 M 9 D 1 M 10 I 1 M 18 D 1 M 5 D 109125 M 16 I 4 M 4 I 1 M 11 I 5 M 8 I 1 M 22 D 2 M 4 D 2 M 34 D 1 M 3 D 1 M 41
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_q_multiple_vulgar(self):
        """Test parsing exn_22_q_multiple_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_q_multiple_vulgar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 6)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 6)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments, check_operations=False)

    def check_vulgar(self, alignments, check_operations=True):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m est2genome comb.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showvulgar yes",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.target.id, "gi|330443482|ref|NC_001134.8|")
        self.assertEqual(alignment.score, 4485)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[560077, 560974], [0, 897]])
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"M"))
        self.assertEqual(
            alignment.format("exonerate", "vulgar"),
            """\
vulgar: gi|296142823|ref|NM_001178508.1| 0 897 + gi|330443482|ref|NC_001134.8| 560077 560974 + 4485 M 897 897
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.target.id, "gi|330443753|ref|NC_001148.4|")
        self.assertEqual(alignment.score, 941)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[492933, 492929, 492929, 492908, 492908, 492857,
                              492857, 492851, 492851, 492840, 492840, 492839,
                              492839, 492813, 492813, 492711, 492711, 492710,
                              492710, 492709, 492709, 492654, 492651, 492646,
                              492646, 492639, 492638, 492634, 492633, 492624,
                              492624, 492621, 492620, 492618, 492616, 492614,
                              492611, 492608, 492607, 492585, 492584, 492575,
                              492575, 492507, 492506, 492503, 492503, 492329,
                              492329, 492321, 492320, 492166, 492163, 492155,
                              492154, 492150, 492148, 492139, 492137, 492136,
                              492135, 492130, 492124, 492071, 492070, 492056,
                              492056, 492033],
                             [     2,      6,      8,     29,     30,     81,
                                  83,     89,     92,    103,    106,    107,
                                 108,    134,    137,    239,    240,    241,
                                 242,    243,    244,    299,    299,    304,
                                 306,    313,    313,    317,    317,    326,
                                 327,    330,    330,    332,    332,    334,
                                 334,    337,    337,    359,    359,    368,
                                 369,    437,    437,    440,    441,    615,
                                 616,    624,    624,    778,    778,    786,
                                 786,    790,    790,    799,    799,    800,
                                 800,    805,    805,    858,    858,    872,
                                 873,    896]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(
                alignment.operations,
                bytearray(
                    b"MIMIMIMIMIMIMIMIMIMIMDMIMDMDMIMDMDMDMDMDMIMDMIMIMDMDMDMDMDMDMDMDMIM"
                ),
            )
        self.assertEqual(
            alignment.format("exonerate", "vulgar"),
            """\
vulgar: gi|296142823|ref|NM_001178508.1| 2 896 + gi|330443753|ref|NC_001148.4| 492933 492033 - 941 M 4 4 G 2 0 M 21 21 G 1 0 M 51 51 G 2 0 M 6 6 G 3 0 M 11 11 G 3 0 M 1 1 G 1 0 M 26 26 G 3 0 M 102 102 G 1 0 M 1 1 G 1 0 M 1 1 G 1 0 M 55 55 G 0 3 M 5 5 G 2 0 M 7 7 G 0 1 M 4 4 G 0 1 M 9 9 G 1 0 M 3 3 G 0 1 M 2 2 G 0 2 M 2 2 G 0 3 M 3 3 G 0 1 M 22 22 G 0 1 M 9 9 G 1 0 M 68 68 G 0 1 M 3 3 G 1 0 M 174 174 G 1 0 M 8 8 G 0 1 M 154 154 G 0 3 M 8 8 G 0 1 M 4 4 G 0 2 M 9 9 G 0 2 M 1 1 G 0 1 M 5 5 G 0 6 M 53 53 G 0 1 M 14 14 G 1 0 M 23 23
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 651)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[267809, 267843, 267843, 267851, 267851, 267852,
                              267852, 267859, 267863, 267867, 267867, 267870,
                              267870, 267880, 267880, 267908, 267908, 267942,
                              267942, 267948, 267950, 267960, 267961, 267975,
                              267976, 267991, 267991, 268003, 268003, 268004,
                              268004, 268007, 268007, 268035, 268035, 268038,
                              268038, 268059, 268059, 268066, 268070, 268072,
                              268074, 268080, 268083, 268096, 268098, 268102,
                              268103, 268114, 268114, 268120, 268120, 268132,
                              268133, 268185, 268185, 268192, 268193, 268248,
                              268249, 268256, 268256, 268269, 268271, 268371,
                              268372, 268378, 268378, 268418, 268418, 268424,
                              268426, 268448, 268450, 300684, 300686, 300698,
                              300699, 300717],
                             [    34,     68,     69,     77,     78,     79,
                                  81,     88,     88,     92,     94,     97,
                                 100,    110,    113,    141,    142,    176,
                                 177,    183,    183,    193,    193,    207,
                                 207,    222,    223,    235,    236,    237,
                                 238,    241,    244,    272,    280,    283,
                                 285,    306,    308,    315,    315,    317,
                                 317,    323,    323,    336,    336,    340,
                                 340,    351,    356,    362,    364,    376,
                                 376,    428,    429,    436,    436,    491,
                                 491,    498,    501,    514,    514,    614,
                                 614,    620,    621,    661,    663,    669,
                                 669,    691,    691,    691,    691,    703,
                                 703,    721]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(
                alignment.operations,
                bytearray(
                    b"MIMIMIMDMIMIMIMIMIMDMDMDMIMIMIMIMIMIMIMDMDMDMDMDMIMIMDMIMDMDMIMDMDMIMIMDM3N5MDM"
                ),
            )
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: gi|296142823|ref|NM_001178508.1| 34 721 + gi|330443520|ref|NC_001136.10| 267809 300717 + 651 M 34 34 G 1 0 M 8 8 G 1 0 M 1 1 G 2 0 M 7 7 G 0 4 M 4 4 G 2 0 M 3 3 G 3 0 M 10 10 G 3 0 M 28 28 G 1 0 M 34 34 G 1 0 M 6 6 G 0 2 M 10 10 G 0 1 M 14 14 G 0 1 M 15 15 G 1 0 M 12 12 G 1 0 M 1 1 G 1 0 M 3 3 G 3 0 M 28 28 G 8 0 M 3 3 G 2 0 M 21 21 G 2 0 M 7 7 G 0 4 M 2 2 G 0 2 M 6 6 G 0 3 M 13 13 G 0 2 M 4 4 G 0 1 M 11 11 G 5 0 M 6 6 G 2 0 M 12 12 G 0 1 M 52 52 G 1 0 M 7 7 G 0 1 M 55 55 G 0 1 M 7 7 G 3 0 M 13 13 G 0 2 M 100 100 G 0 1 M 6 6 G 1 0 M 40 40 G 2 0 M 6 6 G 0 2 M 22 22 3 0 2 I 0 32234 5 0 2 M 12 12 G 0 1 M 18 18
""",
            )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 6150)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"M"))
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: gi|296143771|ref|NM_001180731.1| 0 1230 + gi|330443520|ref|NC_001136.10| 1319275 1318045 - 6150 M 1230 1230
""",
            )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 439)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[ 85010,  85021,  85021,  85036,  85036,  85040,
                               85040,  85041,  85041,  85049,  85049,  85066,
                               85068, 253972, 253974, 253978, 253979, 253987,
                              253987, 253990, 253990, 254023, 254024, 254031,
                              254033, 254135, 254137, 350957, 350959, 350973,
                              350975, 350985, 350985, 350990, 350992, 351002,
                              351002, 351006, 351007, 351027, 351027, 351042,
                              351043, 351048, 351048, 351052, 351054, 473168,
                              473170, 473190, 473195, 473201],
                             [     0,     11,     12,     27,     29,     33,
                                  34,     35,     36,     44,     48,     65,
                                  65,     65,     65,     69,     69,     77,
                                  78,     81,     83,    116,    116,    123,
                                 123,    225,    225,    225,    225,    239,
                                 239,    249,    251,    256,    256,    266,
                                 268,    272,    272,    292,    293,    308,
                                 308,    313,    316,    320,    320,    320,
                                 320,    340,    340,    346]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(
                alignment.operations,
                bytearray(b"MIMIMIMIMIM5N3MDMIMIMDMDM5N3MDMIMDMIMDMIMDMIM5N3MDM"),
            )
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: gi|296143771|ref|NM_001180731.1| 0 346 + gi|330443688|ref|NC_001145.3| 85010 473201 + 439 M 11 11 G 1 0 M 15 15 G 2 0 M 4 4 G 1 0 M 1 1 G 1 0 M 8 8 G 4 0 M 17 17 5 0 2 I 0 168904 3 0 2 M 4 4 G 0 1 M 8 8 G 1 0 M 3 3 G 2 0 M 33 33 G 0 1 M 7 7 G 0 2 M 102 102 5 0 2 I 0 96820 3 0 2 M 14 14 G 0 2 M 10 10 G 2 0 M 5 5 G 0 2 M 10 10 G 2 0 M 4 4 G 0 1 M 20 20 G 1 0 M 15 15 G 0 1 M 5 5 G 3 0 M 4 4 5 0 2 I 0 122114 3 0 2 M 20 20 G 0 5 M 6 6
""",
            )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 263)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[130198, 130184, 130183, 130179, 130179, 130154,
                              130153, 130144, 130138, 130096, 130096, 130080,
                              130078, 130071, 130070, 130067, 130067, 130044,
                              130044, 130038, 130036, 120683, 120681, 120680,
                              120680, 120669, 120668, 120656, 120656, 120647,
                              120646, 120636, 120636, 120618, 120617, 120612,
                              120610,  11489,  11487,  11471,  11471,  11467,
                               11467,  11456,  11456,  11448,  11448,  11426,
                               11424,  11420,  11418,  11384,  11383,  11380,
                               11379,  11338],
                             [    25,     39,     39,     43,     45,     70,
                                  70,     79,     79,    121,    123,    139,
                                 139,    146,    146,    149,    151,    174,
                                 177,    183,    183,    183,    183,    184,
                                 185,    196,    196,    208,    209,    218,
                                 218,    228,    229,    247,    247,    252,
                                 252,    252,    252,    268,    272,    276,
                                 277,    288,    293,    301,    302,    324,
                                 324,    328,    328,    362,    362,    365,
                                 365,    406]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(
                alignment.operations,
                bytearray(b"MDMIMDMDMIMDMDMIMIM3N5MIMDMIMDMIMDM3N5MIMIMIMIMDMDMDMDM"),
            )
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: gi|296143771|ref|NM_001180731.1| 25 406 + gi|330443688|ref|NC_001145.3| 130198 11338 - 263 M 14 14 G 0 1 M 4 4 G 2 0 M 25 25 G 0 1 M 9 9 G 0 6 M 42 42 G 2 0 M 16 16 G 0 2 M 7 7 G 0 1 M 3 3 G 2 0 M 23 23 G 3 0 M 6 6 3 0 2 I 0 9353 5 0 2 M 1 1 G 1 0 M 11 11 G 0 1 M 12 12 G 1 0 M 9 9 G 0 1 M 10 10 G 1 0 M 18 18 G 0 1 M 5 5 3 0 2 I 0 109121 5 0 2 M 16 16 G 4 0 M 4 4 G 1 0 M 11 11 G 5 0 M 8 8 G 1 0 M 22 22 G 0 2 M 4 4 G 0 2 M 34 34 G 0 1 M 3 3 G 0 1 M 41 41
""",
            )
        with self.assertRaises(StopIteration):
            next(alignments)


class Exonerate_coding2coding_fshifts(unittest.TestCase):
    def test_exn_22_m_coding2coding_fshifts_cigar(self):
        """Test parsing exn_22_m_cigar_fshifts.exn)."""
        exn_file = os.path.join("Exonerate", "exn_22_m_cigar_fshifts.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 2)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 2)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)

    def check_cigar(self, alignments):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m coding2coding c2c_frameshift2.fa scer_cad1.fa --showcigar yes --showvulgar no --showalignment no --bestn 3",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.score, 213)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[465, 558, 558, 591, 593, 605, 609, 630],
                             [  0,  93,  94, 127, 127, 139, 139, 160]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 0 160 + gi|296143771|ref|NM_001180731.1| 465 630 + 213 M 93 I 1 M 33 D 2 M 12 D 4 M 21
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.score, 201)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[628, 604, 598, 559, 559, 466],
                             [158, 134, 134,  95,  94,   1]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: gi|296143771|ref|NM_001180731.1| 158 1 - gi|296143771|ref|NM_001180731.1| 628 466 - 201 M 24 D 6 M 39 I 1 M 93
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_coding2coding_fshifts_vulgar(self):
        """Test parsing exn_22_o_vulgar_fshifts.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_o_vulgar_fshifts.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 2)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 2)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments, check_operations=False)

    def check_vulgar(self, alignments, check_operations=True):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m coding2coding c2c_frameshift2.fa scer_cad1.fa --showvulgar yes --showalignment no --bestn 3",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.score, 213)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[465, 558, 558, 591, 593, 605, 609, 630],
                             [  0,  93,  94, 127, 127, 139, 139, 160]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"CFCFCFC"))
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: gi|296143771|ref|NM_001180731.1| 0 160 + gi|296143771|ref|NM_001180731.1| 465 630 + 213 C 93 93 F 1 0 C 33 33 F 0 2 C 12 12 F 0 4 C 21 21
""",
            )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.score, 201)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[628, 604, 598, 559, 559, 466],
                             [158, 134, 134,  95,  94,   1]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"CDCFC"))
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: gi|296143771|ref|NM_001180731.1| 158 1 - gi|296143771|ref|NM_001180731.1| 628 466 - 201 C 24 24 G 0 6 C 39 39 F 1 0 C 93 93
""",
            )
        with self.assertRaises(StopIteration):
            next(alignments)


class Exonerate_protein2dna(unittest.TestCase):
    def test_exn_22_m_protein2dna_cigar(self):
        """Test parsing exonerate output (exn_22_m_protein2dna_cigar.exn)."""
        exn_file = os.path.join("Exonerate", "exn_22_m_protein2dna_cigar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)

    def check_cigar(self, alignments):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m protein2dna ../scer_cad1_prot.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar yes --showvulgar no",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2105)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318048], [0, 409]])
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: sp|P24813|YAP2_YEAST 0 409 . gi|330443520|ref|NC_001136.10| 1319275 1318048 - 2105 M 1227
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 205)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[253991, 254027, 254030, 254270], [28, 40, 40, 120]]),
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: sp|P24813|YAP2_YEAST 28 120 . gi|330443688|ref|NC_001145.3| 253991 254270 + 205 M 36 D 3 M 240
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 116)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[255638, 255743, 255743, 255794], [355, 390, 391, 408]]),
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: sp|P24813|YAP2_YEAST 355 408 . gi|330443688|ref|NC_001145.3| 255638 255794 + 116 M 105 I 1 M 51
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_protein2dna_vulgar(self):
        """Test parsing exonerate output (exn_22_m_protein2dna_vulgar.exn)."""
        exn_file = os.path.join("Exonerate", "exn_22_m_protein2dna_vulgar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments, check_operations=False)

    def check_vulgar(self, alignments, check_operations=True):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m protein2dna ../scer_cad1_prot.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar no --showvulgar yes",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2105)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318048], [0, 409]])
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"M"))
        self.assertEqual(
            alignment.format("exonerate", "vulgar"),
            """\
vulgar: sp|P24813|YAP2_YEAST 0 409 . gi|330443520|ref|NC_001136.10| 1319275 1318048 - 2105 M 409 1227
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 205)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[253991, 254027, 254030, 254270],
                             [    28,     40,     40,    120]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"MDM"))
        self.assertEqual(
            alignment.format("exonerate", "vulgar"),
            """\
vulgar: sp|P24813|YAP2_YEAST 28 120 . gi|330443688|ref|NC_001145.3| 253991 254270 + 205 M 12 36 G 0 3 M 80 240
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 116)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[255638, 255743, 255743, 255794],
                             [   355,    390,    391,    408]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"MIM"))
        self.assertEqual(
            alignment.format("exonerate", "vulgar"),
            """\
vulgar: sp|P24813|YAP2_YEAST 355 408 . gi|330443688|ref|NC_001145.3| 255638 255794 + 116 M 35 105 G 1 0 M 17 51
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)


class Exonerate_protein2dna_fshifts(unittest.TestCase):
    def test_exn_22_m_protein2dna_fshifts_cigar(self):
        """Test parsing exonerate output (exn_22_o_cigar_fshifts2.exn)."""
        exn_file = os.path.join("Exonerate", "exn_22_o_cigar_fshifts2.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 2)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 2)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)

    def check_cigar(self, alignments):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m protein2dna scer_cad1_prot.fa scer_cad1frameshift.fa --showcigar yes --showvulgar no --showalignment no --bestn 3",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.target.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.score, 367)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[216, 345, 347, 455], [330, 373, 373, 409]]),
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: sp|P24813|YAP2_YEAST 330 409 . gi|296143771|ref|NM_001180731.1| 216 455 + 367 M 129 D 2 M 108
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.target.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.score, 322)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[16, 208], [6, 70]]))
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: sp|P24813|YAP2_YEAST 6 70 . gi|296143771|ref|NM_001180731.1| 16 208 + 322 M 192
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_protein2dna_fshifts_vulgar(self):
        """Test parsing exonerate output (exn_22_o_vulgar_fshifts2.exn)."""
        exn_file = os.path.join("Exonerate", "exn_22_o_vulgar_fshifts2.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 2)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 2)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments, check_operations=False)

    def check_vulgar(self, alignments, check_operations=True):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m protein2dna scer_cad1_prot.fa scer_cad1frameshift.fa --showvulgar yes --showalignment no --bestn 3",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.target.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.score, 367)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[216, 345, 347, 455], [330, 373, 373, 409]]),
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"MFM"))
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: sp|P24813|YAP2_YEAST 330 409 . gi|296143771|ref|NM_001180731.1| 216 455 + 367 M 43 129 F 0 2 M 36 108
""",
            )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.target.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.score, 322)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(alignment.coordinates, numpy.array([[16, 208], [6, 70]]))
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"M"))
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: sp|P24813|YAP2_YEAST 6 70 . gi|296143771|ref|NM_001180731.1| 16 208 + 322 M 64 192
""",
            )
        with self.assertRaises(StopIteration):
            next(alignments)


class Exonerate_protein2genome(unittest.TestCase):
    def test_exn_22_m_protein2genome_cigar(self):
        """Test parsing exn_22_m_protein2genome_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_protein2genome_cigar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)

    def check_cigar(self, alignments):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m protein2genome ../scer_cad1_prot.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar yes --showvulgar no",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2105)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318048], [0, 409]])
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: sp|P24813|YAP2_YEAST 0 409 . gi|330443520|ref|NC_001136.10| 1319275 1318048 - 2105 M 1227
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 205)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[253991, 254027, 254030, 254270],
                             [    28,     40,     40,    120]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: sp|P24813|YAP2_YEAST 28 120 . gi|330443688|ref|NC_001145.3| 253991 254270 + 205 M 36 D 3 M 240
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.target.id, "gi|330443590|ref|NC_001140.6|")
        self.assertEqual(alignment.score, 122)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[84646, 84535, 68601, 68450],
                             [   37,    74,    74,   125]])
                # fmt: on
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_protein2genome_vulgar(self):
        """Test parsing exn_22_m_protein2genome_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_protein2genome_vulgar.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 3)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments)
        # For the last alignment, the vulgar string contains a split codon
        # which is not shown in the cigar string.  Writing the alignment with
        # a cigar string and reading it in does not fully regenerate therefore
        # does not regenerate all information provided in the vulgar file.

    def check_vulgar(self, alignments, check_operations=True):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m protein2genome ../scer_cad1_prot.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar no --showvulgar yes",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2105)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318048], [0, 409]])
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"M"))
        self.assertEqual(
            alignment.format("exonerate", "vulgar"),
            """\
vulgar: sp|P24813|YAP2_YEAST 0 409 . gi|330443520|ref|NC_001136.10| 1319275 1318048 - 2105 M 409 1227
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 205)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[253991, 254027, 254030, 254270],
                             [    28,     40,     40,    120]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"MDM"))
        self.assertEqual(
            alignment.format("exonerate", "vulgar"),
            """\
vulgar: sp|P24813|YAP2_YEAST 28 120 . gi|330443688|ref|NC_001145.3| 253991 254270 + 205 M 12 36 G 0 3 M 80 240
""",
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.target.id, "gi|330443590|ref|NC_001140.6|")
        self.assertEqual(alignment.score, 122)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[84646, 84535, 84533, 84531, 68603, 68601,
                              68600, 68450],
                             [   37,    74,    74,    74,    74,    74,
                                 75,   125]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"MS5N3SM"))
        self.assertEqual(
            alignment.format("exonerate", "vulgar"),
            """\
vulgar: sp|P24813|YAP2_YEAST 37 125 . gi|330443590|ref|NC_001140.6| 84646 68450 - 122 M 37 111 S 0 2 5 0 2 I 0 15928 3 0 2 S 1 1 M 50 150
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)


class Exonerate_protein2genome_revcomp_fshifts(unittest.TestCase):
    def test_exn_24_m_protein2genome_revcomp_fshifts_cigar(self):
        """Test parsing exn_24_m_protein2genome_revcomp_fshifts_cigar.exn."""
        exn_file = os.path.join(
            "Exonerate", "exn_24_m_protein2genome_revcomp_fshifts_cigar.exn"
        )
        alignments = Align.parse(exn_file, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)

    def check_cigar(self, alignments):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m protein2genome gene026_baits.fasta gene026_contigs.fasta --showalignment no --showcigar yes --showvulgar no --bestn 2 --refine full",
        )
        self.assertEqual(alignments.metadata["Hostname"], "Michiels-MacBook-Pro.local")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "Morus-gene026")
        self.assertEqual(alignment.target.id, "NODE_2_length_1708_cov_48.590765")
        self.assertEqual(alignment.score, 1308)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[1416, 1380, 1374, 1125, 1125, 1047, 1047,  744,
                               744,  450,  450,  448,  331],
                             [  69,   81,   81,  164,  169,  195,  196,  297,
                               300,  398,  402,  402,  441]])
                # fmt: on
            )
        )
        self.assertEqual(
            alignment.format("exonerate", "cigar"),
            """\
cigar: Morus-gene026 69 441 . NODE_2_length_1708_cov_48.590765 1416 331 - 1308 M 36 D 6 M 249 I 5 M 78 I 1 M 303 I 3 M 294 I 4 D 2 M 117
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_24_m_protein2genome_revcomp_fshifts_vulgar(self):
        """Test parsing exn_24_m_protein2genome_revcomp_fshifts_vulgar.exn."""
        exn_file = os.path.join(
            "Exonerate", "exn_24_m_protein2genome_revcomp_fshifts_vulgar.exn"
        )
        alignments = Align.parse(exn_file, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments, check_operations=False)

    def check_vulgar(self, alignments, check_operations=True):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m protein2genome gene026_baits.fasta gene026_contigs.fasta --showalignment no --showcigar no --showvulgar yes --bestn 2 --refine full",
        )
        self.assertEqual(alignments.metadata["Hostname"], "Michiels-MacBook-Pro.local")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "Morus-gene026")
        self.assertEqual(alignment.target.id, "NODE_2_length_1708_cov_48.590765")
        self.assertEqual(alignment.score, 1308)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[1416, 1380, 1374, 1125, 1125, 1047, 1047,  744,
                               744,  450,  450, 448,  331],
                             [  69,   81,   81,  164,  169,  195,  196,  297,
                               300,  398,  402, 402,  441]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(alignment.operations, bytearray(b"MDMIMIMIMIFM"))
            self.assertEqual(
                alignment.format("exonerate", "vulgar"),
                """\
vulgar: Morus-gene026 69 441 . NODE_2_length_1708_cov_48.590765 1416 331 - 1308 M 12 36 G 0 6 M 83 249 G 5 0 M 26 78 G 1 0 M 101 303 G 3 0 M 98 294 G 4 0 F 0 2 M 39 117
""",
            )
        with self.assertRaises(StopIteration):
            next(alignments)


class Exonerate_protein2genome_met_intron(unittest.TestCase):
    def test_exn_24_protein2genome_met_intron_cigar(self):
        """Test parsing exn_24_m_protein2genome_met_intron_cigar.exn."""
        exn_file = os.path.join(
            "Exonerate", "exn_24_m_protein2genome_met_intron_cigar.exn"
        )
        alignments = Align.parse(exn_file, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "cigar")
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_cigar(alignments)

    def check_cigar(self, alignments):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m protein2genome gene001_baits.fasta gene001_contigs.fasta --showalignment no --showcigar yes --showvulgar no --bestn 1 --refine full",
        )
        self.assertEqual(alignments.metadata["Hostname"], "Michiels-MacBook-Pro.local")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "Morus-gene001")
        self.assertEqual(alignment.target.id, "NODE_1_length_2817_cov_100.387732")
        self.assertEqual(alignment.score, 1978)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[2392, 2281, 2129, 2030, 1921, 1810, 1724, 1421,
                              1198, 1058,  925,  388],
                             [  48,   85,   85,  118,  118,  155,  155,  256,
                               256,  303,  303,  482]])
                # fmt: on
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_24_protein2genome_met_intron_vulgar(self):
        """Test parsing exn_24_m_protein2genome_met_intron_vulgar.exn."""
        exn_file = os.path.join(
            "Exonerate", "exn_24_m_protein2genome_met_intron_vulgar.exn"
        )
        alignments = Align.parse(exn_file, "exonerate")
        self.check_vulgar(alignments)
        alignments = Align.parse(exn_file, "exonerate")
        stream = io.StringIO()
        n = Align.write(alignments, stream, "exonerate", "vulgar")
        self.assertEqual(n, 1)
        stream.seek(0)
        alignments = Align.parse(stream, "exonerate")
        self.check_vulgar(alignments)
        # The vulgar string contains a split codon that is not shown in the
        # cigar string.  Writing the alignment with a cigar string and reading
        # it in does not fully regenerate therefore does not regenerate all
        # information provided in the vulgar file.

    def check_vulgar(self, alignments, check_operations=True):
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m protein2genome gene001_baits.fasta gene001_contigs.fasta --showalignment no --showcigar no --showvulgar yes --bestn 1 --refine full",
        )
        self.assertEqual(alignments.metadata["Hostname"], "Michiels-MacBook-Pro.local")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "Morus-gene001")
        self.assertEqual(alignment.target.id, "NODE_1_length_2817_cov_100.387732")
        self.assertEqual(alignment.score, 1978)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[2392, 2281, 2279, 2131, 2129, 2030, 2028, 1923,
                              1921, 1810, 1808, 1726, 1724, 1421, 1420, 1418,
                              1200, 1198, 1196, 1058, 1056,  927,  925,  388],
                             [  48,   85,   85,   85,   85,  118,  118,  118,
                               118,  155,  155,  155,  155,  256,  256,  256,
                               256,  256,  257,  303,  303,  303,  303,  482]])
                # fmt: on
            )
        )
        if check_operations:
            self.assertEqual(
                alignment.operations, bytearray(b"M5N3M5N3M5N3MS5N3SM5N3M")
            )
        self.assertEqual(
            alignment.format("exonerate", "vulgar"),
            """\
vulgar: Morus-gene001 48 482 . NODE_1_length_2817_cov_100.387732 2392 388 - 1978 M 37 111 5 0 2 I 0 148 3 0 2 M 33 99 5 0 2 I 0 105 3 0 2 M 37 111 5 0 2 I 0 82 3 0 2 M 101 303 S 0 1 5 0 2 I 0 218 3 0 2 S 1 2 M 46 138 5 0 2 I 0 129 3 0 2 M 179 537
""",
        )
        with self.assertRaises(StopIteration):
            next(alignments)


class Exonerate_none(unittest.TestCase):
    def test_exn_22_q_none(self):
        """Test parsing exonerate output (exn_22_q_none.exn)."""
        exn_file = os.path.join("Exonerate", "exn_22_q_none.exn")
        alignments = Align.parse(exn_file, "exonerate")
        self.assertEqual(alignments.metadata["Program"], "exonerate")
        self.assertEqual(
            alignments.metadata["Command line"],
            "exonerate -m est2genome none.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showcigar yes --showvulgar yes",
        )
        self.assertEqual(alignments.metadata["Hostname"], "blackbriar")
        self.assertRaises(StopIteration, next, alignments)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
