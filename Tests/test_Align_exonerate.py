# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for Bio.Align Exonerate parsers."""

import os
import unittest
import warnings

from Bio import BiopythonExperimentalWarning

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    from Bio.Align import exonerate

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.exonerate."
    ) from None


def prettyprint(coordinates):
    n = 0
    for row in coordinates:
        for value in row:
            n = max(n, len(str(value)))
    print()
    fmt = "%%%dd" % n
    for row in coordinates:
        words = []
        for value in row:
            word = fmt % value
            words.append(word)
        line = ", ".join(words)
        print(line)
    

class Exonerate_est2genome(unittest.TestCase):

    def test_exn_22_m_est2genome_cigar(self):
        """Test parsing exn_22_m_est2genome_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_est2genome_cigar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m est2genome ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar yes --showvulgar no")
        self.assertEqual(alignments.hostname, "blackbriar")
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
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_est2genome_vulgar(self):
        """Test parsing exn_22_m_est2genome_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_est2genome_vulgar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m est2genome ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar no --showvulgar yes")
        self.assertEqual(alignments.hostname, "blackbriar")
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
        with self.assertRaises(StopIteration):
            next(alignments)


class Exonerate_affine_local(unittest.TestCase):

    def test_exn_22_m_affine_local_cigar(self):
        """Test parsing exn_22_m_affine_local_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_affine_local_cigar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m affine:local ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar yes --showvulgar no")
        self.assertEqual(alignments.hostname, "Michiels-MacBook-Pro.local")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 6150)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
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
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_affine_local_vulgar(self):
        """Test parsing exn_22_m_affine_local_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_affine_local_vulgar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m affine:local ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar no --showvulgar yes")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 6150)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
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
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443715|ref|NC_001146.8|")
        self.assertEqual(alignment.score, 219)
        with self.assertRaises(StopIteration):
            next(alignments)


class Exonerate_cdna2genome(unittest.TestCase):

    def test_exn_22_m_cdna2genome_cigar(self):
        """Test parsing exn_22_m_cdna2genome_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_cdna2genome_cigar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m cdna2genome ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar yes --showvulgar no")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 6146)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(

                alignment.coordinates,
                numpy.array([[1319275, 1319274, 1319271, 1318045],
                             [      0,       1,       4,    1230]]
                )
            )
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
                numpy.array([[1318045, 1318174, 1318177, 1319275],
                             [   1230,    1101,    1098,       0]])
            )
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
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_cdna2genome_vulgar(self):
        """Test parsing exn_22_m_cdna2genome_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_cdna2genome_vulgar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m cdna2genome ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar no --showvulgar yes")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 6146)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[1319275, 1319274, 1319271, 1318045],
                             [      0,       1,       4,    1230]])
            )
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
                numpy.array([[1318045, 1318174, 1318177, 1319275],
                             [   1230,    1101,    1098,       0]])
            )
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
        with self.assertRaises(StopIteration):
            next(alignments)


class Exonerate_coding2coding(unittest.TestCase):

    def test_exn_22_m_coding2coding_cigar(self):
        """Test parsing exn_22_m_coding2coding_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_coding2coding_cigar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m coding2coding ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar yes --showvulgar no")
        self.assertEqual(alignments.hostname, "blackbriar")
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
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 116)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[255638, 255743, 255743, 255794],
                             [  1065,   1170,   1173,   1224]])
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_coding2coding_vulgar(self):
        """Test parsing exn_22_m_coding2coding_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_coding2coding_vulgar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m coding2coding ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar no --showvulgar yes")
        self.assertEqual(alignments.hostname, "blackbriar")
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
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 116)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[255638, 255743, 255743, 255794],
                             [  1065,   1170,   1173,   1224]])
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)


class Exonerate_coding2genome(unittest.TestCase):

    def test_exn_22_m_coding2genome_cigar(self):
        """Test parsing exn_22_m_coding2genome_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_coding2genome_cigar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m coding2genome ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar yes --showvulgar no")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2151)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1318047, 1319274], [1228, 1]])
                # fmt: on
            )
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
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 116)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[255638, 255743, 255743, 255794],
                             [  1065,   1170,   1173,   1224]])
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_coding2genome_vulgar(self):
        """Test parsing exn_22_m_coding2genome_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_coding2genome_vulgar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m coding2genome ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showcigar no --showvulgar yes")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2151)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1318047, 1319274], [1228, 1]])
                # fmt: on
            )
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
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.score, 116)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[255638, 255743, 255743, 255794],
                             [  1065,   1170,   1173,   1224]])
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)


class Exonerate_dna2protein(unittest.TestCase):

    def test_exn_22_m_dna2protein_cigar(self):
        """Test parsing exn_22_m_dna2protein_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_dna2protein_cigar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate --showcigar yes --showvulgar no --showalignment no nuc2.fa pro.fa")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "dna")
        self.assertEqual(alignment.target.id, "protein")
        self.assertEqual(alignment.score, 105)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[313, 344], [0, 93]])
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_dna2protein_vulgar(self):
        """Test parsing exn_22_m_dna2protein_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_dna2protein_vulgar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate --showcigar no --showvulgar yes --showalignment no nuc2.fa pro.fa")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "dna")
        self.assertEqual(alignment.target.id, "protein")
        self.assertEqual(alignment.score, 105)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[313, 344], [0, 93]])
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)


class Exonerate_genome2genome(unittest.TestCase):

    def test_exn_22_m_genome2genome_cigar(self):
        """Test parsing exn_22_o_vulgar_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_o_vulgar_cigar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m genome2genome ../intron.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showvulgar yes --showcigar yes")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sacCer3_dna")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2641)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[1319997, 1319971, 1319968, 1319468],
                             [    529,     503,     500,       0]])
            )
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
                numpy.array([[1319997, 1319971, 1319968, 1319468],
                             [    529,     503,     500,       0]])
            )
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
                 numpy.array([[1319468, 1319558, 1319561, 1319997],
                              [      0,      90,      93,     529]])
            )
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
                 numpy.array([[1319468, 1319558, 1319561, 1319997],
                              [      0,      90,      93,     529]])
            )
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
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_genome2genome_vulgar(self):
        """Test parsing exn_22_o_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_o_vulgar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m genome2genome ../intron.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showalignment no --showvulgar yes --showcigar no")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sacCer3_dna")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.score, 2641)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[1319997, 1319971, 1319968, 1319468],
                             [    529,     503,     500,       0]])
            )
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
                 numpy.array([[1319468, 1319558, 1319561, 1319997],
                              [      0,      90,      93,     529]])
            )
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
        with self.assertRaises(StopIteration):
            next(alignments)


class ExonerateTextCases(unittest.TestCase):

    def test_exn_22_m_ungapped_cigar(self):
        """Test parsing exn_22_m_ungapped_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_ungapped_cigar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m ungapped ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showcigar no --showvulgar no")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.annotations["model"], "ungapped:dna2dna")
        self.assertEqual(alignment.score, 6150)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(repr(alignment.target.seq), "Seq({1318045: 'CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTAT...CAT'}, length=1319275)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        self.assertEqual(alignment[0], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(alignment[1], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.annotations["model"], "ungapped:dna2dna")
        self.assertEqual(alignment.score, 233)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({121: 'TTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAG...AGA'}, length=236)")
        self.assertEqual(repr(alignment.target.seq), "Seq({254031: 'TGGATCCTGAAACTAAGCAGAAGAGGACTGCCCAAAATCGGGCCGCTCAAAGAG...AGA'}, length=254146)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[254031, 254146], [242, 357]])
            )
        )
        self.assertEqual(alignment[0], "TGGATCCTGAAACTAAGCAGAAGAGGACTGCCCAAAATCGGGCCGCTCAAAGAGCTTTTAGGGAACGTAAGGAGAGGAAGATGAAGGAATTGGAGAAGAAGGTACAAAGTTTAGA")
        self.assertEqual(alignment[1], "TTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGA")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence")
        self.assertEqual(alignment.annotations["model"], "ungapped:dna2dna")
        self.assertEqual(alignment.score, 151)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({1098: 'CCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAG...AGA'}, length=1166)")
        self.assertEqual(repr(alignment.target.seq), "Seq({255671: 'CCGAAATACTCAGATATTGATGTCGATGGTTTATGTTCCGAGCTAATGGCAAAG...AGA'}, length=255739)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[255671, 255739], [2196, 2264]])
            )
        )
        self.assertEqual(alignment[0], "CCGAAATACTCAGATATTGATGTCGATGGTTTATGTTCCGAGCTAATGGCAAAGGCAAAATGTTCAGA")
        self.assertEqual(alignment[1], "CCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGA")
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_ungapped_vulgar(self):
        """Test parsing exn_22_m_ungapped_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_ungapped_vulgar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m ungapped ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showcigar no --showvulgar no")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "ungapped:dna2dna")
        self.assertEqual(alignment.score, 6150)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(repr(alignment.target.seq), "Seq({1318045: 'CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTAT...CAT'}, length=1319275)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        self.assertEqual(alignment[0], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(alignment[1], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence")
        self.assertEqual(alignment.annotations["model"], "ungapped:dna2dna")
        self.assertEqual(alignment.score, 233)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({121: 'TTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAG...AGA'}, length=236)")
        self.assertEqual(repr(alignment.target.seq), "Seq({254031: 'TGGATCCTGAAACTAAGCAGAAGAGGACTGCCCAAAATCGGGCCGCTCAAAGAG...AGA'}, length=254146)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[254031, 254146], [242, 357]])
            )
        )
        self.assertEqual(alignment[0], "TGGATCCTGAAACTAAGCAGAAGAGGACTGCCCAAAATCGGGCCGCTCAAAGAGCTTTTAGGGAACGTAAGGAGAGGAAGATGAAGGAATTGGAGAAGAAGGTACAAAGTTTAGA")
        self.assertEqual(alignment[1], "TTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGA")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence")
        self.assertEqual(alignment.annotations["model"], "ungapped:dna2dna")
        self.assertEqual(alignment.score, 151)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({1098: 'CCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAG...AGA'}, length=1166)")
        self.assertEqual(repr(alignment.target.seq), "Seq({255671: 'CCGAAATACTCAGATATTGATGTCGATGGTTTATGTTCCGAGCTAATGGCAAAG...AGA'}, length=255739)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[255671, 255739], [2196, 2264]])
            )
        )
        self.assertEqual(alignment[0], "CCGAAATACTCAGATATTGATGTCGATGGTTTATGTTCCGAGCTAATGGCAAAGGCAAAATGTTCAGA")
        self.assertEqual(alignment[1], "CCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGA")
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_ungapped_trans_cigar(self):
        """Test parsing exn_22_m_ungapped_trans_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_ungapped_trans_cigar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m ungapped:trans ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showcigar no --showvulgar no")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds:[revcomp]")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence")
        self.assertEqual(alignment.annotations["model"], "ungapped:codon")
        self.assertEqual(alignment.score, 2151)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({1: 'TGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGA...TGT'}, length=1228)")
        self.assertEqual(repr(alignment.target.seq), "Seq({1318047: 'ACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTATTT...CCA'}, length=1319274)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1318047, 1319274], [1228, 1]])
            )
        )
        self.assertEqual(alignment[0], "ACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTATTTTGCAGTCATCTGTACATTTTGCCTTGATTATTAATTCGCTGCATAAATCATCTATGTCCAACGATGAATATTTTGGTAGGGAGGAGATCTCTTCGAGAATGTGGTAGCAAGACGCTTTTCCATTCGATATGGCGCTGCTAATGTGTTCAATTGCTTCAGATTGGGTTCGAATTGTTTTAGTATGAGAAGCAGCCACCGGAAGCACCGAAGCTAGAGATTTCTGATTCAACAGGTGGCTTGTTATTAAGTCGCATTTGATCTCCTTATTTGAACATGGTTTGGTCCCGATACACTTGCGGTTGCAAATATTATTACACTGGCTCCCTACAGCAAGGGGAAGCACGCTGGGGAATAAAGTATCATCTTCTGTCTCGTATGGAAGATACTTTTGGTCAATTTCACATTGACATGCTGACGTATCTATTTTTTCAGCGGTCAGGCGTTCAGAGTGCTGAGGAAAGTAGTTATGAATGTCGTAATCTAAAGCTTCTAGACCAGTCTGTCTATCATTGGAGACAAGATCGAGACTGAAATCTGGCATTTCACTACCATTAATACTTTCGTCACCGATGGTAACGGCACCTGTTCGATTTCCGGGGGGAGAATCCCACGAATTGGTTAAACTGGGAGCGTTCCATGAGCCTATGTTCCACATAACTTTCTCATTCTCTTGAGTGTGCGGCTGTACTTCTATACTTGTGTTCACAGTCATGTTTTCATCTGAATTAGGCGATGGCAATTCCTTTGCAGCCTTGCTAACTGCTGTACTTGTTCCTTTTTCGTTTTCCCTTTTCTGTTGTTCTTGCAGATCATCGAGGAAGGCTAATATTCTTTCATCATCAGAATTCTTAGCTCTATATTTTGTAATTTCCGACAGTAAACTTTTTAAAGAACATAGTAAAAAGTCCGTGGTAGTCTTATTCTGCGCATCTTTCTGTTCTAGTAACTCTACCCTCTCTTGCAAACTCTTCATTTTGGCTTCTTTCCTATCTCGGAACGCTCGTTGAGCTGCCCTATTCTGGGCAGTCCTCCTACTCTTAGCTTCAGAGTCAATCCTTTTTCTGCCTGGTCTGCCCACCTTCCTCTTAGGTGTTCCATCTTTATTTAGCAACATTTGCTTCTTCATGTCACCTGCTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCA")
        self.assertEqual(alignment[1], "ACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTATTTTGCAGTCATCTGTACATTTTGCCTTGATTATTAATTCGCTGCATAAATCATCTATGTCCAACGATGAATATTTTGGTAGGGAGGAGATCTCTTCGAGAATGTGGTAGCAAGACGCTTTTCCATTCGATATGGCGCTGCTAATGTGTTCAATTGCTTCAGATTGGGTTCGAATTGTTTTAGTATGAGAAGCAGCCACCGGAAGCACCGAAGCTAGAGATTTCTGATTCAACAGGTGGCTTGTTATTAAGTCGCATTTGATCTCCTTATTTGAACATGGTTTGGTCCCGATACACTTGCGGTTGCAAATATTATTACACTGGCTCCCTACAGCAAGGGGAAGCACGCTGGGGAATAAAGTATCATCTTCTGTCTCGTATGGAAGATACTTTTGGTCAATTTCACATTGACATGCTGACGTATCTATTTTTTCAGCGGTCAGGCGTTCAGAGTGCTGAGGAAAGTAGTTATGAATGTCGTAATCTAAAGCTTCTAGACCAGTCTGTCTATCATTGGAGACAAGATCGAGACTGAAATCTGGCATTTCACTACCATTAATACTTTCGTCACCGATGGTAACGGCACCTGTTCGATTTCCGGGGGGAGAATCCCACGAATTGGTTAAACTGGGAGCGTTCCATGAGCCTATGTTCCACATAACTTTCTCATTCTCTTGAGTGTGCGGCTGTACTTCTATACTTGTGTTCACAGTCATGTTTTCATCTGAATTAGGCGATGGCAATTCCTTTGCAGCCTTGCTAACTGCTGTACTTGTTCCTTTTTCGTTTTCCCTTTTCTGTTGTTCTTGCAGATCATCGAGGAAGGCTAATATTCTTTCATCATCAGAATTCTTAGCTCTATATTTTGTAATTTCCGACAGTAAACTTTTTAAAGAACATAGTAAAAAGTCCGTGGTAGTCTTATTCTGCGCATCTTTCTGTTCTAGTAACTCTACCCTCTCTTGCAAACTCTTCATTTTGGCTTCTTTCCTATCTCGGAACGCTCGTTGAGCTGCCCTATTCTGGGCAGTCCTCCTACTCTTAGCTTCAGAGTCAATCCTTTTTCTGCCTGGTCTGCCCACCTTCCTCTTAGGTGTTCCATCTTTATTTAGCAACATTTGCTTCTTCATGTCACCTGCTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCA")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "ungapped:codon")
        self.assertEqual(alignment.score, 2106)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(repr(alignment.target.seq), "Seq({1318045: 'CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTAT...CAT'}, length=1319275)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        self.assertEqual(alignment[0], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(alignment[1], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds:[revcomp]")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence")
        self.assertEqual(alignment.annotations["model"], "ungapped:codon")
        self.assertEqual(alignment.score, 2072)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(repr(alignment.target.seq), "Seq({1318045: 'CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTAT...CAT'}, length=1319275)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1318045, 1319275], [1230, 0]])
            )
        )
        self.assertEqual(alignment[0], "CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTATTTTGCAGTCATCTGTACATTTTGCCTTGATTATTAATTCGCTGCATAAATCATCTATGTCCAACGATGAATATTTTGGTAGGGAGGAGATCTCTTCGAGAATGTGGTAGCAAGACGCTTTTCCATTCGATATGGCGCTGCTAATGTGTTCAATTGCTTCAGATTGGGTTCGAATTGTTTTAGTATGAGAAGCAGCCACCGGAAGCACCGAAGCTAGAGATTTCTGATTCAACAGGTGGCTTGTTATTAAGTCGCATTTGATCTCCTTATTTGAACATGGTTTGGTCCCGATACACTTGCGGTTGCAAATATTATTACACTGGCTCCCTACAGCAAGGGGAAGCACGCTGGGGAATAAAGTATCATCTTCTGTCTCGTATGGAAGATACTTTTGGTCAATTTCACATTGACATGCTGACGTATCTATTTTTTCAGCGGTCAGGCGTTCAGAGTGCTGAGGAAAGTAGTTATGAATGTCGTAATCTAAAGCTTCTAGACCAGTCTGTCTATCATTGGAGACAAGATCGAGACTGAAATCTGGCATTTCACTACCATTAATACTTTCGTCACCGATGGTAACGGCACCTGTTCGATTTCCGGGGGGAGAATCCCACGAATTGGTTAAACTGGGAGCGTTCCATGAGCCTATGTTCCACATAACTTTCTCATTCTCTTGAGTGTGCGGCTGTACTTCTATACTTGTGTTCACAGTCATGTTTTCATCTGAATTAGGCGATGGCAATTCCTTTGCAGCCTTGCTAACTGCTGTACTTGTTCCTTTTTCGTTTTCCCTTTTCTGTTGTTCTTGCAGATCATCGAGGAAGGCTAATATTCTTTCATCATCAGAATTCTTAGCTCTATATTTTGTAATTTCCGACAGTAAACTTTTTAAAGAACATAGTAAAAAGTCCGTGGTAGTCTTATTCTGCGCATCTTTCTGTTCTAGTAACTCTACCCTCTCTTGCAAACTCTTCATTTTGGCTTCTTTCCTATCTCGGAACGCTCGTTGAGCTGCCCTATTCTGGGCAGTCCTCCTACTCTTAGCTTCAGAGTCAATCCTTTTTCTGCCTGGTCTGCCCACCTTCCTCTTAGGTGTTCCATCTTTATTTAGCAACATTTGCTTCTTCATGTCACCTGCTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCAT")
        self.assertEqual(alignment[1], "CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTATTTTGCAGTCATCTGTACATTTTGCCTTGATTATTAATTCGCTGCATAAATCATCTATGTCCAACGATGAATATTTTGGTAGGGAGGAGATCTCTTCGAGAATGTGGTAGCAAGACGCTTTTCCATTCGATATGGCGCTGCTAATGTGTTCAATTGCTTCAGATTGGGTTCGAATTGTTTTAGTATGAGAAGCAGCCACCGGAAGCACCGAAGCTAGAGATTTCTGATTCAACAGGTGGCTTGTTATTAAGTCGCATTTGATCTCCTTATTTGAACATGGTTTGGTCCCGATACACTTGCGGTTGCAAATATTATTACACTGGCTCCCTACAGCAAGGGGAAGCACGCTGGGGAATAAAGTATCATCTTCTGTCTCGTATGGAAGATACTTTTGGTCAATTTCACATTGACATGCTGACGTATCTATTTTTTCAGCGGTCAGGCGTTCAGAGTGCTGAGGAAAGTAGTTATGAATGTCGTAATCTAAAGCTTCTAGACCAGTCTGTCTATCATTGGAGACAAGATCGAGACTGAAATCTGGCATTTCACTACCATTAATACTTTCGTCACCGATGGTAACGGCACCTGTTCGATTTCCGGGGGGAGAATCCCACGAATTGGTTAAACTGGGAGCGTTCCATGAGCCTATGTTCCACATAACTTTCTCATTCTCTTGAGTGTGCGGCTGTACTTCTATACTTGTGTTCACAGTCATGTTTTCATCTGAATTAGGCGATGGCAATTCCTTTGCAGCCTTGCTAACTGCTGTACTTGTTCCTTTTTCGTTTTCCCTTTTCTGTTGTTCTTGCAGATCATCGAGGAAGGCTAATATTCTTTCATCATCAGAATTCTTAGCTCTATATTTTGTAATTTCCGACAGTAAACTTTTTAAAGAACATAGTAAAAAGTCCGTGGTAGTCTTATTCTGCGCATCTTTCTGTTCTAGTAACTCTACCCTCTCTTGCAAACTCTTCATTTTGGCTTCTTTCCTATCTCGGAACGCTCGTTGAGCTGCCCTATTCTGGGCAGTCCTCCTACTCTTAGCTTCAGAGTCAATCCTTTTTCTGCCTGGTCTGCCCACCTTCCTCTTAGGTGTTCCATCTTTATTTAGCAACATTTGCTTCTTCATGTCACCTGCTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCAT")
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_ungapped_trans_vulgar(self):
        """Test parsing exn_22_m_ungapped_trans_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_ungapped_trans_vulgar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m ungapped:trans ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showcigar no --showvulgar no")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds:[revcomp]")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence")
        self.assertEqual(alignment.annotations["model"], "ungapped:codon")
        self.assertEqual(alignment.score, 2151)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({1: 'TGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGA...TGT'}, length=1228)")
        self.assertEqual(repr(alignment.target.seq), "Seq({1318047: 'ACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTATTT...CCA'}, length=1319274)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1318047, 1319274], [1228, 1]])
            )
        )
        self.assertEqual(alignment[0], "ACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTATTTTGCAGTCATCTGTACATTTTGCCTTGATTATTAATTCGCTGCATAAATCATCTATGTCCAACGATGAATATTTTGGTAGGGAGGAGATCTCTTCGAGAATGTGGTAGCAAGACGCTTTTCCATTCGATATGGCGCTGCTAATGTGTTCAATTGCTTCAGATTGGGTTCGAATTGTTTTAGTATGAGAAGCAGCCACCGGAAGCACCGAAGCTAGAGATTTCTGATTCAACAGGTGGCTTGTTATTAAGTCGCATTTGATCTCCTTATTTGAACATGGTTTGGTCCCGATACACTTGCGGTTGCAAATATTATTACACTGGCTCCCTACAGCAAGGGGAAGCACGCTGGGGAATAAAGTATCATCTTCTGTCTCGTATGGAAGATACTTTTGGTCAATTTCACATTGACATGCTGACGTATCTATTTTTTCAGCGGTCAGGCGTTCAGAGTGCTGAGGAAAGTAGTTATGAATGTCGTAATCTAAAGCTTCTAGACCAGTCTGTCTATCATTGGAGACAAGATCGAGACTGAAATCTGGCATTTCACTACCATTAATACTTTCGTCACCGATGGTAACGGCACCTGTTCGATTTCCGGGGGGAGAATCCCACGAATTGGTTAAACTGGGAGCGTTCCATGAGCCTATGTTCCACATAACTTTCTCATTCTCTTGAGTGTGCGGCTGTACTTCTATACTTGTGTTCACAGTCATGTTTTCATCTGAATTAGGCGATGGCAATTCCTTTGCAGCCTTGCTAACTGCTGTACTTGTTCCTTTTTCGTTTTCCCTTTTCTGTTGTTCTTGCAGATCATCGAGGAAGGCTAATATTCTTTCATCATCAGAATTCTTAGCTCTATATTTTGTAATTTCCGACAGTAAACTTTTTAAAGAACATAGTAAAAAGTCCGTGGTAGTCTTATTCTGCGCATCTTTCTGTTCTAGTAACTCTACCCTCTCTTGCAAACTCTTCATTTTGGCTTCTTTCCTATCTCGGAACGCTCGTTGAGCTGCCCTATTCTGGGCAGTCCTCCTACTCTTAGCTTCAGAGTCAATCCTTTTTCTGCCTGGTCTGCCCACCTTCCTCTTAGGTGTTCCATCTTTATTTAGCAACATTTGCTTCTTCATGTCACCTGCTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCA")
        self.assertEqual(alignment[1], "ACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTATTTTGCAGTCATCTGTACATTTTGCCTTGATTATTAATTCGCTGCATAAATCATCTATGTCCAACGATGAATATTTTGGTAGGGAGGAGATCTCTTCGAGAATGTGGTAGCAAGACGCTTTTCCATTCGATATGGCGCTGCTAATGTGTTCAATTGCTTCAGATTGGGTTCGAATTGTTTTAGTATGAGAAGCAGCCACCGGAAGCACCGAAGCTAGAGATTTCTGATTCAACAGGTGGCTTGTTATTAAGTCGCATTTGATCTCCTTATTTGAACATGGTTTGGTCCCGATACACTTGCGGTTGCAAATATTATTACACTGGCTCCCTACAGCAAGGGGAAGCACGCTGGGGAATAAAGTATCATCTTCTGTCTCGTATGGAAGATACTTTTGGTCAATTTCACATTGACATGCTGACGTATCTATTTTTTCAGCGGTCAGGCGTTCAGAGTGCTGAGGAAAGTAGTTATGAATGTCGTAATCTAAAGCTTCTAGACCAGTCTGTCTATCATTGGAGACAAGATCGAGACTGAAATCTGGCATTTCACTACCATTAATACTTTCGTCACCGATGGTAACGGCACCTGTTCGATTTCCGGGGGGAGAATCCCACGAATTGGTTAAACTGGGAGCGTTCCATGAGCCTATGTTCCACATAACTTTCTCATTCTCTTGAGTGTGCGGCTGTACTTCTATACTTGTGTTCACAGTCATGTTTTCATCTGAATTAGGCGATGGCAATTCCTTTGCAGCCTTGCTAACTGCTGTACTTGTTCCTTTTTCGTTTTCCCTTTTCTGTTGTTCTTGCAGATCATCGAGGAAGGCTAATATTCTTTCATCATCAGAATTCTTAGCTCTATATTTTGTAATTTCCGACAGTAAACTTTTTAAAGAACATAGTAAAAAGTCCGTGGTAGTCTTATTCTGCGCATCTTTCTGTTCTAGTAACTCTACCCTCTCTTGCAAACTCTTCATTTTGGCTTCTTTCCTATCTCGGAACGCTCGTTGAGCTGCCCTATTCTGGGCAGTCCTCCTACTCTTAGCTTCAGAGTCAATCCTTTTTCTGCCTGGTCTGCCCACCTTCCTCTTAGGTGTTCCATCTTTATTTAGCAACATTTGCTTCTTCATGTCACCTGCTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCA")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "ungapped:codon")
        self.assertEqual(alignment.score, 2106)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(repr(alignment.target.seq), "Seq({1318045: 'CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTAT...CAT'}, length=1319275)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        self.assertEqual(alignment[0], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(alignment[1], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds:[revcomp]")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence")
        self.assertEqual(alignment.annotations["model"], "ungapped:codon")
        self.assertEqual(alignment.score, 2072)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(repr(alignment.target.seq), "Seq({1318045: 'CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTAT...CAT'}, length=1319275)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1318045, 1319275], [1230, 0]])
            )
        )
        self.assertEqual(alignment[0], "CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTATTTTGCAGTCATCTGTACATTTTGCCTTGATTATTAATTCGCTGCATAAATCATCTATGTCCAACGATGAATATTTTGGTAGGGAGGAGATCTCTTCGAGAATGTGGTAGCAAGACGCTTTTCCATTCGATATGGCGCTGCTAATGTGTTCAATTGCTTCAGATTGGGTTCGAATTGTTTTAGTATGAGAAGCAGCCACCGGAAGCACCGAAGCTAGAGATTTCTGATTCAACAGGTGGCTTGTTATTAAGTCGCATTTGATCTCCTTATTTGAACATGGTTTGGTCCCGATACACTTGCGGTTGCAAATATTATTACACTGGCTCCCTACAGCAAGGGGAAGCACGCTGGGGAATAAAGTATCATCTTCTGTCTCGTATGGAAGATACTTTTGGTCAATTTCACATTGACATGCTGACGTATCTATTTTTTCAGCGGTCAGGCGTTCAGAGTGCTGAGGAAAGTAGTTATGAATGTCGTAATCTAAAGCTTCTAGACCAGTCTGTCTATCATTGGAGACAAGATCGAGACTGAAATCTGGCATTTCACTACCATTAATACTTTCGTCACCGATGGTAACGGCACCTGTTCGATTTCCGGGGGGAGAATCCCACGAATTGGTTAAACTGGGAGCGTTCCATGAGCCTATGTTCCACATAACTTTCTCATTCTCTTGAGTGTGCGGCTGTACTTCTATACTTGTGTTCACAGTCATGTTTTCATCTGAATTAGGCGATGGCAATTCCTTTGCAGCCTTGCTAACTGCTGTACTTGTTCCTTTTTCGTTTTCCCTTTTCTGTTGTTCTTGCAGATCATCGAGGAAGGCTAATATTCTTTCATCATCAGAATTCTTAGCTCTATATTTTGTAATTTCCGACAGTAAACTTTTTAAAGAACATAGTAAAAAGTCCGTGGTAGTCTTATTCTGCGCATCTTTCTGTTCTAGTAACTCTACCCTCTCTTGCAAACTCTTCATTTTGGCTTCTTTCCTATCTCGGAACGCTCGTTGAGCTGCCCTATTCTGGGCAGTCCTCCTACTCTTAGCTTCAGAGTCAATCCTTTTTCTGCCTGGTCTGCCCACCTTCCTCTTAGGTGTTCCATCTTTATTTAGCAACATTTGCTTCTTCATGTCACCTGCTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCAT")
        self.assertEqual(alignment[1], "CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTATTTTGCAGTCATCTGTACATTTTGCCTTGATTATTAATTCGCTGCATAAATCATCTATGTCCAACGATGAATATTTTGGTAGGGAGGAGATCTCTTCGAGAATGTGGTAGCAAGACGCTTTTCCATTCGATATGGCGCTGCTAATGTGTTCAATTGCTTCAGATTGGGTTCGAATTGTTTTAGTATGAGAAGCAGCCACCGGAAGCACCGAAGCTAGAGATTTCTGATTCAACAGGTGGCTTGTTATTAAGTCGCATTTGATCTCCTTATTTGAACATGGTTTGGTCCCGATACACTTGCGGTTGCAAATATTATTACACTGGCTCCCTACAGCAAGGGGAAGCACGCTGGGGAATAAAGTATCATCTTCTGTCTCGTATGGAAGATACTTTTGGTCAATTTCACATTGACATGCTGACGTATCTATTTTTTCAGCGGTCAGGCGTTCAGAGTGCTGAGGAAAGTAGTTATGAATGTCGTAATCTAAAGCTTCTAGACCAGTCTGTCTATCATTGGAGACAAGATCGAGACTGAAATCTGGCATTTCACTACCATTAATACTTTCGTCACCGATGGTAACGGCACCTGTTCGATTTCCGGGGGGAGAATCCCACGAATTGGTTAAACTGGGAGCGTTCCATGAGCCTATGTTCCACATAACTTTCTCATTCTCTTGAGTGTGCGGCTGTACTTCTATACTTGTGTTCACAGTCATGTTTTCATCTGAATTAGGCGATGGCAATTCCTTTGCAGCCTTGCTAACTGCTGTACTTGTTCCTTTTTCGTTTTCCCTTTTCTGTTGTTCTTGCAGATCATCGAGGAAGGCTAATATTCTTTCATCATCAGAATTCTTAGCTCTATATTTTGTAATTTCCGACAGTAAACTTTTTAAAGAACATAGTAAAAAGTCCGTGGTAGTCTTATTCTGCGCATCTTTCTGTTCTAGTAACTCTACCCTCTCTTGCAAACTCTTCATTTTGGCTTCTTTCCTATCTCGGAACGCTCGTTGAGCTGCCCTATTCTGGGCAGTCCTCCTACTCTTAGCTTCAGAGTCAATCCTTTTTCTGCCTGGTCTGCCCACCTTCCTCTTAGGTGTTCCATCTTTATTTAGCAACATTTGCTTCTTCATGTCACCTGCTAAATATATTTGCTGACCTTTCCGAAGGATATTGCCCAT")
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_ner_cigar(self):
        """Test parsing exonerate output (exn_22_m_ner_cigar.exn)."""
        exn_file = os.path.join("Exonerate", "exn_22_m_ner_cigar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m ner ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showcigar no --showvulgar no")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "NER:affine:local:dna2dna")
        self.assertEqual(alignment.score, 6150)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(repr(alignment.target.seq), "Seq({1318045: 'CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTAT...CAT'}, length=1319275)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        self.assertEqual(alignment[0], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(alignment[1], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443681|ref|NC_001144.5|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XII, complete sequence")
        self.assertEqual(alignment.annotations["model"], "NER:affine:local:dna2dna")
        self.assertEqual(alignment.score, 502)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({110: 'CAGAAAA', 148: 'CTGCCCAGAAT', 169: 'AACGAGCGTTCCG', 184: 'ATAGGAAAGAAGC', 227: 'GTTACTAGAACAGAAAG', 255: 'AAGACT', 266: 'GGACTTTTTACTATGTTC', 296: 'ACTGTCGGAAAT', 309: 'ACAAAATATAGAGCTAAGAATTCTGATGAT', 708: 'AACTACTT', 728: 'TGAACG', 744: 'GAAAAAATAGATACGTCAGCATG', 778: 'TTGACCAAAAGTATCTTCCATACGAG', 816: 'ACTTT', 832: 'TGCTTCCCCTTGC', 889: 'GGACCAA', 915: 'ATCAAATGCGACTTAA', 940: 'ACCTGT', 953: 'GAAATC', 960: 'CTAGCTTC', 979: 'TGGCTGCTTCTCAT', 1006: 'GAACCCAAT', 1016: 'TGAAGC', 1024: 'TTGAACACATTAGCAG', 1047: 'TCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAA', 1092: 'TCCCT', 1099: 'CAAAATATTCATCGTTGGA', 1121: 'AGATGATTT', 1133: 'CAGCGA', 1156: 'AATGTACAG', 1167: 'GACTGCAAAA', 1187: 'AGCTCGCGACTTA', 1210: 'TGGTTAGACAGCTCCTGTAG'}, length=1230)")
        self.assertEqual(repr(alignment.target.seq), "Seq({297910: 'CTGAAAA', 297946: 'CCGCCCAAAGT', 297970: 'AACTGGAGTTCCG', 297992: 'ATTGGAAGATGC', 298019: 'GTTATCAAGAACAATAAAG', 298049: 'AAGACT', 298065: 'GCACTATTCACTGTCCTC', 298095: 'ACATTCTCAAAT', 298117: 'ACAAAAAAACAGCTGATTGATTATTTGATGAA', 318437: 'AGCTTCTT', 318454: 'TGAAAG', 318476: 'GAATAAATACATGGCACCTTG', 318503: 'TTGACATTCAATGCCTTCTAAAGAG', 318565: 'ACTTT', 318611: 'TTCCTCCCCTTCC', 318630: 'GGACGAA', 318639: 'ATCAAATCAACCAAA', 318659: 'ACCTGT', 318681: 'GAAATC', 318691: 'CTAGTTTC', 318708: 'TGTCTCATCTCAT', 318725: 'GAATCCAGT', 318742: 'TGTAGC', 318755: 'TTGACCACATTCACCAG', 318787: 'TCGAATGCAAAGAAGCTAGCTGAACATTATCGAA', 318834: 'TCCCT', 318851: 'CAAAATTTTCTTCATAGAA', 318889: 'AGATATTT', 318913: 'CAGCTA', 318924: 'ACTGTACCG', 318944: 'GACTTCAAAA', 318959: 'AGCTATAGACTCA', 318974: 'TGATAGGACAGCTCCTGTAG'}, length=318994)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[297910, 297917, 297917, 297946, 297957, 297957,
                              297970, 297983, 297983, 297992, 297997, 297997,
                              298004, 298004, 298019, 298023, 298024, 298031,
                              298032, 298038, 298038, 298049, 298055, 298055,
                              298065, 298083, 298083, 298095, 298107, 298107,
                              298117, 298123, 298123, 298133, 298135, 298139,
                              298140, 298149, 298149, 318437, 318445, 318445,
                              318454, 318460, 318460, 318476, 318488, 318488,
                              318497, 318497, 318503, 318523, 318523, 318528,
                              318528, 318565, 318570, 318570, 318611, 318624,
                              318624, 318630, 318637, 318637, 318639, 318646,
                              318646, 318654, 318654, 318659, 318665, 318665,
                              318681, 318687, 318687, 318691, 318699, 318699,
                              318708, 318713, 318713, 318721, 318721, 318725,
                              318734, 318734, 318742, 318748, 318748, 318755,
                              318766, 318767, 318772, 318772, 318787, 318793,
                              318793, 318808, 318808, 318821, 318821, 318834,
                              318839, 318839, 318851, 318870, 318870, 318889,
                              318893, 318893, 318897, 318897, 318913, 318919,
                              318919, 318924, 318933, 318933, 318944, 318954,
                              318954, 318959, 318972, 318972, 318974, 318994],
                             [   110,    117,    148,    148,    159,    169,
                                 169,    182,    184,    184,    189,    190,
                                 197,    227,    227,    231,    231,    238,
                                 238,    244,    255,    255,    261,    266,
                                 266,    284,    296,    296,    308,    309,
                                 309,    315,    316,    326,    326,    330,
                                 330,    339,    708,    708,    716,    728,
                                 728,    734,    744,    744,    756,    758,
                                 767,    778,    778,    798,    799,    804,
                                 816,    816,    821,    832,    832,    845,
                                 889,    889,    896,    915,    915,    922,
                                 923,    931,    940,    940,    946,    953,
                                 953,    959,    960,    960,    968,    979,
                                 979,    984,    985,    993,   1006,   1006,
                                1015,   1016,   1016,   1022,   1024,   1024,
                                1035,   1035,   1040,   1047,   1047,   1053,
                                1054,   1069,   1070,   1083,   1092,   1092,
                                1097,   1099,   1099,   1118,   1121,   1121,
                                1125,   1126,   1130,   1133,   1133,   1139,
                                1156,   1156,   1165,   1167,   1167,   1177,
                                1187,   1187,   1200,   1210,   1210,   1230]])
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence")
        self.assertEqual(alignment.annotations["model"], "NER:affine:local:dna2dna")
        self.assertEqual(alignment.score, 440)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({509: 'TGAGA', 537: 'TGGAACGCTC', 567: 'GATTCTCCCCCCGGAAATCGAACAGGTG', 607: 'GTGACGAAAG', 636: 'CCAGATTTCAGTCT', 667: 'ATAGACAGACTGGTC', 683: 'AGAAGCTTTAGATTACGACATTCATAACTACTTTC', 737: 'GACCGCTGA', 748: 'AAATAGATACGTCAGCAT', 777: 'ATTGACCAAAAGTA', 797: 'ATACGAGACAGA', 815: 'TACTTTATTCCCCAGCGT', 843: 'GCTGT', 861: 'AATAATATTTGCAACCGCAAG', 902: 'TTCAAATAAGGAGATCAA', 927: 'TTAATAACAAG', 957: 'TCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATAC', 1010: 'CCAATCTGAA', 1022: 'AATTGAACA', 1050: 'AATGGAAAAGCGTCTTGCTACCACATTCT', 1099: 'CAAAATATT', 1122: 'GATGATTTATGCAGCGAA', 1147: 'TCAAGGC', 1167: 'GACTGCAAAATAGTAGTCAAAGCTC'}, length=1192)")
        self.assertEqual(repr(alignment.target.seq), "Seq({183946: 'TGAGA', 183977: 'TGAAACTGTC', 184002: 'GAAGCTCTTTCAGAAAAGCCAACTCGTG', 184044: 'GTGAAGACAG', 184066: 'CCGGATTACAATCT', 184092: 'AGAGAAAGGCAAGTC', 184111: 'AGAAACTGGCAGATCTTGGACTTGAAAAGGAGTTTC', 184169: 'GGCCGCTGA', 184190: 'AAATGCAGATCATCAGCAT', 184211: 'ATTACGAAAACTA', 184239: 'AAACCAGAAAGA', 184271: 'TACTAGGTTACCAGCTT', 184299: 'GCTGT', 184315: 'AAGAAGACATGGAAACCCAGG', 184346: 'TTCAAAGAAGAAGAAGAA', 184366: 'TTAAAAAGAAG', 184393: 'TGTCTACTAGAGAACTTGCCGGCAAGGTTGGGTCTATAC', 184439: 'CAAATCTGGA', 184463: 'AATTGGAAA', 184485: 'ATTGGCAAAGGTGCCGAAACCACATTTT', 184515: 'CAAGATGTT', 184542: 'GATGATGCATCCTCAGCGGA', 184565: 'TCTAGGC', 184578: 'GACGGTAAAATAGTAGTCACACCTC'}, length=184603)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[183946, 183951, 183951, 183977, 183987, 183987,
                              184002, 184030, 184030, 184044, 184054, 184054,
                              184066, 184080, 184080, 184092, 184107, 184107,
                              184111, 184118, 184119, 184125, 184126, 184132,
                              184132, 184147, 184147, 184169, 184178, 184178,
                              184190, 184194, 184196, 184200, 184200, 184209,
                              184209, 184211, 184214, 184214, 184224, 184224,
                              184239, 184251, 184251, 184271, 184280, 184280,
                              184288, 184288, 184299, 184304, 184304, 184315,
                              184336, 184336, 184346, 184364, 184364, 184366,
                              184377, 184377, 184393, 184399, 184399, 184410,
                              184411, 184415, 184417, 184428, 184428, 184432,
                              184432, 184439, 184449, 184449, 184463, 184472,
                              184472, 184485, 184495, 184495, 184513, 184513,
                              184515, 184524, 184524, 184542, 184552, 184554,
                              184562, 184562, 184565, 184572, 184572, 184578,
                              184603],
                             [   509,    514,    537,    537,    547,    567,
                                 567,    595,    607,    607,    617,    636,
                                 636,    650,    667,    667,    682,    683,
                                 683,    690,    690,    696,    696,    702,
                                 703,    718,    737,    737,    746,    748,
                                 748,    752,    752,    756,    757,    766,
                                 777,    777,    780,    781,    791,    797,
                                 797,    809,    815,    815,    824,    825,
                                 833,    843,    843,    848,    861,    861,
                                 882,    902,    902,    920,    927,    927,
                                 938,    957,    957,    963,    964,    975,
                                 975,    979,    979,    990,    991,    995,
                                1010,   1010,   1020,   1022,   1022,   1031,
                                1050,   1050,   1060,   1061,   1079,   1099,
                                1099,   1108,   1122,   1122,   1132,   1132,
                                1140,   1147,   1147,   1154,   1167,   1167,
                                1192]])
                # fmt: on
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_ner_vulgar(self):
        """Test parsing exonerate output (exn_22_m_ner_vulgar.exn)."""
        exn_file = os.path.join("Exonerate", "exn_22_m_ner_vulgar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m ner ../scer_cad1.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showcigar no --showvulgar no")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "NER:affine:local:dna2dna")
        self.assertEqual(alignment.score, 6150)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(repr(alignment.target.seq), "Seq({1318045: 'CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTAT...CAT'}, length=1319275)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        self.assertEqual(alignment[0], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(alignment[1], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443681|ref|NC_001144.5|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XII, complete sequence")
        self.assertEqual(alignment.annotations["model"], "NER:affine:local:dna2dna")
        self.assertEqual(alignment.score, 502)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({110: 'CAGAAAA', 148: 'CTGCCCAGAAT', 169: 'AACGAGCGTTCCG', 184: 'ATAGGAAAGAAGC', 227: 'GTTACTAGAACAGAAAG', 255: 'AAGACT', 266: 'GGACTTTTTACTATGTTC', 296: 'ACTGTCGGAAAT', 309: 'ACAAAATATAGAGCTAAGAATTCTGATGAT', 708: 'AACTACTT', 728: 'TGAACG', 744: 'GAAAAAATAGATACGTCAGCATG', 778: 'TTGACCAAAAGTATCTTCCATACGAG', 816: 'ACTTT', 832: 'TGCTTCCCCTTGC', 889: 'GGACCAA', 915: 'ATCAAATGCGACTTAA', 940: 'ACCTGT', 953: 'GAAATC', 960: 'CTAGCTTC', 979: 'TGGCTGCTTCTCAT', 1006: 'GAACCCAAT', 1016: 'TGAAGC', 1024: 'TTGAACACATTAGCAG', 1047: 'TCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAA', 1092: 'TCCCT', 1099: 'CAAAATATTCATCGTTGGA', 1121: 'AGATGATTT', 1133: 'CAGCGA', 1156: 'AATGTACAG', 1167: 'GACTGCAAAA', 1187: 'AGCTCGCGACTTA', 1210: 'TGGTTAGACAGCTCCTGTAG'}, length=1230)")
        self.assertEqual(repr(alignment.target.seq), "Seq({297910: 'CTGAAAA', 297946: 'CCGCCCAAAGT', 297970: 'AACTGGAGTTCCG', 297992: 'ATTGGAAGATGC', 298019: 'GTTATCAAGAACAATAAAG', 298049: 'AAGACT', 298065: 'GCACTATTCACTGTCCTC', 298095: 'ACATTCTCAAAT', 298117: 'ACAAAAAAACAGCTGATTGATTATTTGATGAA', 318437: 'AGCTTCTT', 318454: 'TGAAAG', 318476: 'GAATAAATACATGGCACCTTG', 318503: 'TTGACATTCAATGCCTTCTAAAGAG', 318565: 'ACTTT', 318611: 'TTCCTCCCCTTCC', 318630: 'GGACGAA', 318639: 'ATCAAATCAACCAAA', 318659: 'ACCTGT', 318681: 'GAAATC', 318691: 'CTAGTTTC', 318708: 'TGTCTCATCTCAT', 318725: 'GAATCCAGT', 318742: 'TGTAGC', 318755: 'TTGACCACATTCACCAG', 318787: 'TCGAATGCAAAGAAGCTAGCTGAACATTATCGAA', 318834: 'TCCCT', 318851: 'CAAAATTTTCTTCATAGAA', 318889: 'AGATATTT', 318913: 'CAGCTA', 318924: 'ACTGTACCG', 318944: 'GACTTCAAAA', 318959: 'AGCTATAGACTCA', 318974: 'TGATAGGACAGCTCCTGTAG'}, length=318994)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[297910, 297917, 297917, 297946, 297957, 297957,
                              297970, 297983, 297983, 297992, 297997, 297997,
                              298004, 298004, 298019, 298023, 298024, 298031,
                              298032, 298038, 298038, 298049, 298055, 298055,
                              298065, 298083, 298083, 298095, 298107, 298107,
                              298117, 298123, 298123, 298133, 298135, 298139,
                              298140, 298149, 298149, 318437, 318445, 318445,
                              318454, 318460, 318460, 318476, 318488, 318488,
                              318497, 318497, 318503, 318523, 318523, 318528,
                              318528, 318565, 318570, 318570, 318611, 318624,
                              318624, 318630, 318637, 318637, 318639, 318646,
                              318646, 318654, 318654, 318659, 318665, 318665,
                              318681, 318687, 318687, 318691, 318699, 318699,
                              318708, 318713, 318713, 318721, 318721, 318725,
                              318734, 318734, 318742, 318748, 318748, 318755,
                              318766, 318767, 318772, 318772, 318787, 318793,
                              318793, 318808, 318808, 318821, 318821, 318834,
                              318839, 318839, 318851, 318870, 318870, 318889,
                              318893, 318893, 318897, 318897, 318913, 318919,
                              318919, 318924, 318933, 318933, 318944, 318954,
                              318954, 318959, 318972, 318972, 318974, 318994],
                             [   110,    117,    148,    148,    159,    169,
                                 169,    182,    184,    184,    189,    190,
                                 197,    227,    227,    231,    231,    238,
                                 238,    244,    255,    255,    261,    266,
                                 266,    284,    296,    296,    308,    309,
                                 309,    315,    316,    326,    326,    330,
                                 330,    339,    708,    708,    716,    728,
                                 728,    734,    744,    744,    756,    758,
                                 767,    778,    778,    798,    799,    804,
                                 816,    816,    821,    832,    832,    845,
                                 889,    889,    896,    915,    915,    922,
                                 923,    931,    940,    940,    946,    953,
                                 953,    959,    960,    960,    968,    979,
                                 979,    984,    985,    993,   1006,   1006,
                                1015,   1016,   1016,   1022,   1024,   1024,
                                1035,   1035,   1040,   1047,   1047,   1053,
                                1054,   1069,   1070,   1083,   1092,   1092,
                                1097,   1099,   1099,   1118,   1121,   1121,
                                1125,   1126,   1130,   1133,   1133,   1139,
                                1156,   1156,   1165,   1167,   1167,   1177,
                                1187,   1187,   1200,   1210,   1210,   1230]])
                # fmt: on
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence")
        self.assertEqual(alignment.annotations["model"], "NER:affine:local:dna2dna")
        self.assertEqual(alignment.score, 440)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({509: 'TGAGA', 537: 'TGGAACGCTC', 567: 'GATTCTCCCCCCGGAAATCGAACAGGTG', 607: 'GTGACGAAAG', 636: 'CCAGATTTCAGTCT', 667: 'ATAGACAGACTGGTC', 683: 'AGAAGCTTTAGATTACGACATTCATAACTACTTTC', 737: 'GACCGCTGA', 748: 'AAATAGATACGTCAGCAT', 777: 'ATTGACCAAAAGTA', 797: 'ATACGAGACAGA', 815: 'TACTTTATTCCCCAGCGT', 843: 'GCTGT', 861: 'AATAATATTTGCAACCGCAAG', 902: 'TTCAAATAAGGAGATCAA', 927: 'TTAATAACAAG', 957: 'TCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATAC', 1010: 'CCAATCTGAA', 1022: 'AATTGAACA', 1050: 'AATGGAAAAGCGTCTTGCTACCACATTCT', 1099: 'CAAAATATT', 1122: 'GATGATTTATGCAGCGAA', 1147: 'TCAAGGC', 1167: 'GACTGCAAAATAGTAGTCAAAGCTC'}, length=1192)")
        self.assertEqual(repr(alignment.target.seq), "Seq({183946: 'TGAGA', 183977: 'TGAAACTGTC', 184002: 'GAAGCTCTTTCAGAAAAGCCAACTCGTG', 184044: 'GTGAAGACAG', 184066: 'CCGGATTACAATCT', 184092: 'AGAGAAAGGCAAGTC', 184111: 'AGAAACTGGCAGATCTTGGACTTGAAAAGGAGTTTC', 184169: 'GGCCGCTGA', 184190: 'AAATGCAGATCATCAGCAT', 184211: 'ATTACGAAAACTA', 184239: 'AAACCAGAAAGA', 184271: 'TACTAGGTTACCAGCTT', 184299: 'GCTGT', 184315: 'AAGAAGACATGGAAACCCAGG', 184346: 'TTCAAAGAAGAAGAAGAA', 184366: 'TTAAAAAGAAG', 184393: 'TGTCTACTAGAGAACTTGCCGGCAAGGTTGGGTCTATAC', 184439: 'CAAATCTGGA', 184463: 'AATTGGAAA', 184485: 'ATTGGCAAAGGTGCCGAAACCACATTTT', 184515: 'CAAGATGTT', 184542: 'GATGATGCATCCTCAGCGGA', 184565: 'TCTAGGC', 184578: 'GACGGTAAAATAGTAGTCACACCTC'}, length=184603)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                # fmt: off
# flake8: noqa
                numpy.array([[183946, 183951, 183951, 183977, 183987, 183987,
                              184002, 184030, 184030, 184044, 184054, 184054,
                              184066, 184080, 184080, 184092, 184107, 184107,
                              184111, 184118, 184119, 184125, 184126, 184132,
                              184132, 184147, 184147, 184169, 184178, 184178,
                              184190, 184194, 184196, 184200, 184200, 184209,
                              184209, 184211, 184214, 184214, 184224, 184224,
                              184239, 184251, 184251, 184271, 184280, 184280,
                              184288, 184288, 184299, 184304, 184304, 184315,
                              184336, 184336, 184346, 184364, 184364, 184366,
                              184377, 184377, 184393, 184399, 184399, 184410,
                              184411, 184415, 184417, 184428, 184428, 184432,
                              184432, 184439, 184449, 184449, 184463, 184472,
                              184472, 184485, 184495, 184495, 184513, 184513,
                              184515, 184524, 184524, 184542, 184552, 184554,
                              184562, 184562, 184565, 184572, 184572, 184578,
                              184603],
                             [   509,    514,    537,    537,    547,    567,
                                 567,    595,    607,    607,    617,    636,
                                 636,    650,    667,    667,    682,    683,
                                 683,    690,    690,    696,    696,    702,
                                 703,    718,    737,    737,    746,    748,
                                 748,    752,    752,    756,    757,    766,
                                 777,    777,    780,    781,    791,    797,
                                 797,    809,    815,    815,    824,    825,
                                 833,    843,    843,    848,    861,    861,
                                 882,    902,    902,    920,    927,    927,
                                 938,    957,    957,    963,    964,    975,
                                 975,    979,    979,    990,    991,    995,
                                1010,   1010,   1020,   1022,   1022,   1031,
                                1050,   1050,   1060,   1061,   1079,   1099,
                                1099,   1108,   1122,   1122,   1132,   1132,
                                1140,   1147,   1147,   1154,   1167,   1167,
                                1192]])
                # fmt: on
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_q_multiple_cigar(self):
        """Test parsing exn_22_q_multiple_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_q_multiple_cigar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m est2genome comb.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showcigar yes --showvulgar yes")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cdc28p (CDC28) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443482|ref|NC_001134.8|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome II, complete sequence")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 4485)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGAGGGTATTCCAAAGGACCAACCGTTAGGAGCTGATATTGTTAAGAAGTTTATGATGCAACTTTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAAAGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGGAAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAGTGGTACCAAGTCTAGATCCACGCGGTATTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCCAGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATAA")
        self.assertEqual(repr(alignment.target.seq), "Seq({560077: 'ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACA...TAA'}, length=560974)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[560077, 560974], [0, 897]])
            )
        )
        self.assertEqual(alignment[0], "ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGAGGGTATTCCAAAGGACCAACCGTTAGGAGCTGATATTGTTAAGAAGTTTATGATGCAACTTTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAAAGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGGAAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAGTGGTACCAAGTCTAGATCCACGCGGTATTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCCAGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATAA")
        self.assertEqual(alignment[1], "ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGAGGGTATTCCAAAGGACCAACCGTTAGGAGCTGATATTGTTAAGAAGTTTATGATGCAACTTTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAAAGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGGAAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAGTGGTACCAAGTCTAGATCCACGCGGTATTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCCAGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATAA")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cdc28p (CDC28) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443482|ref|NC_001134.8|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome II, complete sequence")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 4485)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGAGGGTATTCCAAAGGACCAACCGTTAGGAGCTGATATTGTTAAGAAGTTTATGATGCAACTTTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAAAGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGGAAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAGTGGTACCAAGTCTAGATCCACGCGGTATTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCCAGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATAA")
        self.assertEqual(repr(alignment.target.seq), "Seq({560077: 'ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACA...TAA'}, length=560974)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[560077, 560974], [0, 897]])
            )
        )
        self.assertEqual(alignment[0], "ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGAGGGTATTCCAAAGGACCAACCGTTAGGAGCTGATATTGTTAAGAAGTTTATGATGCAACTTTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAAAGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGGAAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAGTGGTACCAAGTCTAGATCCACGCGGTATTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCCAGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATAA")
        self.assertEqual(alignment[1], "ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGAGGGTATTCCAAAGGACCAACCGTTAGGAGCTGATATTGTTAAGAAGTTTATGATGCAACTTTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAAAGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGGAAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAGTGGTACCAAGTCTAGATCCACGCGGTATTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCCAGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATAA")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cdc28p (CDC28) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443482|ref|NC_001134.8|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome II, complete sequence")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 4485)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGAGGGTATTCCAAAGGACCAACCGTTAGGAGCTGATATTGTTAAGAAGTTTATGATGCAACTTTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAAAGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGGAAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAGTGGTACCAAGTCTAGATCCACGCGGTATTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCCAGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATAA")
        self.assertEqual(repr(alignment.target.seq), "Seq({560077: 'ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACA...TAA'}, length=560974)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[560077, 560974], [0, 897]])
            )
        )
        self.assertEqual(alignment[0], "ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGAGGGTATTCCAAAGGACCAACCGTTAGGAGCTGATATTGTTAAGAAGTTTATGATGCAACTTTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAAAGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGGAAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAGTGGTACCAAGTCTAGATCCACGCGGTATTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCCAGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATAA")
        self.assertEqual(alignment[1], "ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGAGGGTATTCCAAAGGACCAACCGTTAGGAGCTGATATTGTTAAGAAGTTTATGATGCAACTTTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAAAGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGGAAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAGTGGTACCAAGTCTAGATCCACGCGGTATTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCCAGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATAA")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cdc28p (CDC28) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443753|ref|NC_001148.4|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XVI, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 941)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({2: 'GAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATA...ATA'}, length=896)")
        self.assertEqual(repr(alignment.target.seq), "Seq({492033: 'TAGTACTCTGCAAACCAAGGGTGATGCAGAGCCTGCTTGGCGCTCAGCCTCATA...CTC'}, length=492933)")
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
        self.assertEqual(alignment[0], "GAGC--TCTATGAACAGATTTAAGCAG-TTAGAAAAGCTTGGCAATGGTACGTATGCCACAGTGTACAAGGGACTGAAC--AAAACC---ACAGGGGTATA---T-GTTGCCCTGAAAGAGGTAAAACTGGA---TTCAGAGGAAGGTACACCCTCTACGGCCATCCGTGAGATCTCCCTAATGAAAGAATTGAAACATGAGAACATTGTTAGACTTTATGACGTTATTCACACAGA-G-A-ACAAGTTGACTTTGGTTTTTGAATTCATGGACAACGATTTAAAGAAATACATGGATTCCCGCA--CCGTGGGCAACACACCAAGAGG-GCTAGAACTAAACTTGGTTAAATACTTCCAGTGGCAACTACTGCAAGGG-CTGGCCTTTTGCCATGAAAACAAGATTCTCCACCGTGATTTAAAACCTCAAAACCTATTAATCAACAAGAGA-GGCCAGTTGAAATTGGGTGATTTCGGTCTGGCCCGTGCTTTCGGTATTCCGGTCAACACATTTTCAAGCGAAGTCGTAACGTTGTGGTACCGTGCTCCTGATGTGCTAATGGGTTCTAGGACGTACTCCACATCCATTGATATATGGTCGTGTGGGTGCATTCTTGCGGAAATG-ATAACGGGTAAGCCTTTGTTTCCTGGCACCAACGACGAAGAACAACTGAAATTGATCTTCGACATCATGGGCACTCCTAATGAGTCCCTATGGCCCAGTGTAACAAAGTTACCCAAATACAACCCAAATATCCAGCAACGACCACCAAGAGACCTACGTCAAGTATTGCAACCACACACCAAAGAACCGCTAGACGGGAATCTCATGGATTTCTTACACGGACTCTTGCAACTTAATCCGGATATGAGGCTGAGCGCCAAGCAGGCTCTGCAT-CACCCTTGGTTTGCAGAGTACTA")
        self.assertEqual(alignment[1], "GAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGA---GGGTATTCCAAAGG-ACCA-ACCGTTAGGAGCT-GA--TA---TTG-TTAAGAAGTTTATGATGCAACT-TTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAA-AGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGG-AAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAG---TGGTACCA-AGTC--TAGATCCAC--G-CGGTA------TTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCC-AGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATA")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cdc28p (CDC28) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443753|ref|NC_001148.4|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XVI, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 941)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({2: 'GAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATA...ATA'}, length=896)")
        self.assertEqual(repr(alignment.target.seq), "Seq({492033: 'TAGTACTCTGCAAACCAAGGGTGATGCAGAGCCTGCTTGGCGCTCAGCCTCATA...CTC'}, length=492933)")
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
        self.assertEqual(alignment[0], "GAGC--TCTATGAACAGATTTAAGCAG-TTAGAAAAGCTTGGCAATGGTACGTATGCCACAGTGTACAAGGGACTGAAC--AAAACC---ACAGGGGTATA---T-GTTGCCCTGAAAGAGGTAAAACTGGA---TTCAGAGGAAGGTACACCCTCTACGGCCATCCGTGAGATCTCCCTAATGAAAGAATTGAAACATGAGAACATTGTTAGACTTTATGACGTTATTCACACAGA-G-A-ACAAGTTGACTTTGGTTTTTGAATTCATGGACAACGATTTAAAGAAATACATGGATTCCCGCA--CCGTGGGCAACACACCAAGAGG-GCTAGAACTAAACTTGGTTAAATACTTCCAGTGGCAACTACTGCAAGGG-CTGGCCTTTTGCCATGAAAACAAGATTCTCCACCGTGATTTAAAACCTCAAAACCTATTAATCAACAAGAGA-GGCCAGTTGAAATTGGGTGATTTCGGTCTGGCCCGTGCTTTCGGTATTCCGGTCAACACATTTTCAAGCGAAGTCGTAACGTTGTGGTACCGTGCTCCTGATGTGCTAATGGGTTCTAGGACGTACTCCACATCCATTGATATATGGTCGTGTGGGTGCATTCTTGCGGAAATG-ATAACGGGTAAGCCTTTGTTTCCTGGCACCAACGACGAAGAACAACTGAAATTGATCTTCGACATCATGGGCACTCCTAATGAGTCCCTATGGCCCAGTGTAACAAAGTTACCCAAATACAACCCAAATATCCAGCAACGACCACCAAGAGACCTACGTCAAGTATTGCAACCACACACCAAAGAACCGCTAGACGGGAATCTCATGGATTTCTTACACGGACTCTTGCAACTTAATCCGGATATGAGGCTGAGCGCCAAGCAGGCTCTGCAT-CACCCTTGGTTTGCAGAGTACTA")
        self.assertEqual(alignment[1], "GAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGA---GGGTATTCCAAAGG-ACCA-ACCGTTAGGAGCT-GA--TA---TTG-TTAAGAAGTTTATGATGCAACT-TTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAA-AGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGG-AAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAG---TGGTACCA-AGTC--TAGATCCAC--G-CGGTA------TTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCC-AGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATA")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cdc28p (CDC28) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443753|ref|NC_001148.4|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XVI, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 941)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({2: 'GAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATA...ATA'}, length=896)")
        self.assertEqual(repr(alignment.target.seq), "Seq({492033: 'TAGTACTCTGCAAACCAAGGGTGATGCAGAGCCTGCTTGGCGCTCAGCCTCATA...CTC'}, length=492933)")
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
        self.assertEqual(alignment[0], "GAGC--TCTATGAACAGATTTAAGCAG-TTAGAAAAGCTTGGCAATGGTACGTATGCCACAGTGTACAAGGGACTGAAC--AAAACC---ACAGGGGTATA---T-GTTGCCCTGAAAGAGGTAAAACTGGA---TTCAGAGGAAGGTACACCCTCTACGGCCATCCGTGAGATCTCCCTAATGAAAGAATTGAAACATGAGAACATTGTTAGACTTTATGACGTTATTCACACAGA-G-A-ACAAGTTGACTTTGGTTTTTGAATTCATGGACAACGATTTAAAGAAATACATGGATTCCCGCA--CCGTGGGCAACACACCAAGAGG-GCTAGAACTAAACTTGGTTAAATACTTCCAGTGGCAACTACTGCAAGGG-CTGGCCTTTTGCCATGAAAACAAGATTCTCCACCGTGATTTAAAACCTCAAAACCTATTAATCAACAAGAGA-GGCCAGTTGAAATTGGGTGATTTCGGTCTGGCCCGTGCTTTCGGTATTCCGGTCAACACATTTTCAAGCGAAGTCGTAACGTTGTGGTACCGTGCTCCTGATGTGCTAATGGGTTCTAGGACGTACTCCACATCCATTGATATATGGTCGTGTGGGTGCATTCTTGCGGAAATG-ATAACGGGTAAGCCTTTGTTTCCTGGCACCAACGACGAAGAACAACTGAAATTGATCTTCGACATCATGGGCACTCCTAATGAGTCCCTATGGCCCAGTGTAACAAAGTTACCCAAATACAACCCAAATATCCAGCAACGACCACCAAGAGACCTACGTCAAGTATTGCAACCACACACCAAAGAACCGCTAGACGGGAATCTCATGGATTTCTTACACGGACTCTTGCAACTTAATCCGGATATGAGGCTGAGCGCCAAGCAGGCTCTGCAT-CACCCTTGGTTTGCAGAGTACTA")
        self.assertEqual(alignment[1], "GAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGA---GGGTATTCCAAAGG-ACCA-ACCGTTAGGAGCT-GA--TA---TTG-TTAAGAAGTTTATGATGCAACT-TTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAA-AGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGG-AAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAG---TGGTACCA-AGTC--TAGATCCAC--G-CGGTA------TTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCC-AGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATA")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cdc28p (CDC28) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 651)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({34: 'AGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGAC...TCT'}, length=721)")
        self.assertEqual(repr(alignment.target.seq), "Seq({267809: 'AGAAAGTTGGTGAGGGTACTTATGCGGTTGTTTACTTGGGTTGTCAACACTCTA...CCT', 300684: 'AACAAACCGAGATAGAATGGCCAGATATTCTCT'}, length=300717)")
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
        self.assertEqual(alignment.query.annotations["splicesites"], (745,))
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cdc28p (CDC28) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 651)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({34: 'AGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGAC...TCT'}, length=721)")
        self.assertEqual(repr(alignment.target.seq), "Seq({267809: 'AGAAAGTTGGTGAGGGTACTTATGCGGTTGTTTACTTGGGTTGTCAACACTCTA...CCT', 300684: 'AACAAACCGAGATAGAATGGCCAGATATTCTCT'}, length=300717)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
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
                # fmt: off
# flake8: noqa
            )
        )
        self.assertEqual(alignment.query.annotations["splicesites"], (745,))
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cdc28p (CDC28) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 651)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({34: 'AGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGAC...TCT'}, length=721)")
        self.assertEqual(repr(alignment.target.seq), "Seq({267809: 'AGAAAGTTGGTGAGGGTACTTATGCGGTTGTTTACTTGGGTTGTCAACACTCTA...CCT', 300684: 'AACAAACCGAGATAGAATGGCCAGATATTCTCT'}, length=300717)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
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
                # fmt: off
# flake8: noqa
            )
        )
        self.assertEqual(alignment.query.annotations["splicesites"], (745,))
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 6150)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(repr(alignment.target.seq), "Seq({1318045: 'CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTAT...CAT'}, length=1319275)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        self.assertEqual(alignment[0], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(alignment[1], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 6150)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(repr(alignment.target.seq), "Seq({1318045: 'CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTAT...CAT'}, length=1319275)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        self.assertEqual(alignment[0], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(alignment[1], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 6150)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(repr(alignment.target.seq), "Seq({1318045: 'CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTAT...CAT'}, length=1319275)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        self.assertEqual(alignment[0], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(alignment[1], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 439)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAA")
        self.assertEqual(repr(alignment.target.seq), "Seq({85010: 'ATGGTGAACCTCTTCAAGACGGTCAGAATAATCAACAGGATGAAGAAGCAAAAGATGT', 253972: 'TGGCGAGGATAGCGAGCAACCGAAGAAGAAGGGTAGCAAAACTAGCAAAAAGCA...ACA', 350957: 'TGGGGTGATTATATCATTTCTGGATGAGGAATACCTGAAGACCACTCTTCATTA...GTT', 473168: 'AGATGGAAGAATTCTGATAATGCTGTAAAAGAA'}, length=473201)")
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
        self.assertEqual(alignment.query.annotations["splicesites"], (448, 448, 448))
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 439)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAA")
        self.assertEqual(repr(alignment.target.seq), "Seq({85010: 'ATGGTGAACCTCTTCAAGACGGTCAGAATAATCAACAGGATGAAGAAGCAAAAGATGT', 253972: 'TGGCGAGGATAGCGAGCAACCGAAGAAGAAGGGTAGCAAAACTAGCAAAAAGCA...ACA', 350957: 'TGGGGTGATTATATCATTTCTGGATGAGGAATACCTGAAGACCACTCTTCATTA...GTT', 473168: 'AGATGGAAGAATTCTGATAATGCTGTAAAAGAA'}, length=473201)")
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
        self.assertEqual(alignment.query.annotations["splicesites"], (448, 448, 448))
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 439)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAA")
        self.assertEqual(repr(alignment.target.seq), "Seq({85010: 'ATGGTGAACCTCTTCAAGACGGTCAGAATAATCAACAGGATGAAGAAGCAAAAGATGT', 253972: 'TGGCGAGGATAGCGAGCAACCGAAGAAGAAGGGTAGCAAAACTAGCAAAAAGCA...ACA', 350957: 'TGGGGTGATTATATCATTTCTGGATGAGGAATACCTGAAGACCACTCTTCATTA...GTT', 473168: 'AGATGGAAGAATTCTGATAATGCTGTAAAAGAA'}, length=473201)")
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
        self.assertEqual(alignment.query.annotations["splicesites"], (448, 448, 448))
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 263)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({25: 'AGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATG...CAA'}, length=406)")
        self.assertEqual(repr(alignment.target.seq), "Seq({11338: 'TTGTATCTACTTCGTTTTGGCTCTTATAATTTTTTAGTGGCCTAAATTTTCAAA...TGT', 120610: 'TTCTTCTTCTTCTTCTTGTTCTTCTTCTCCATCATTATCTTCTTCATCTTCTTC...CGA', 130036: 'AGTCCGAATCTTTTTGCTTGGCTCCTCTGGGGTCACTTCAACCTTATTAATTTC...GGT'}, length=130198)")
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
        self.assertEqual(alignment.query.annotations["splicesites"], (459, 459))
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 263)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({25: 'AGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATG...CAA'}, length=406)")
        self.assertEqual(repr(alignment.target.seq), "Seq({11338: 'TTGTATCTACTTCGTTTTGGCTCTTATAATTTTTTAGTGGCCTAAATTTTCAAA...TGT', 120610: 'TTCTTCTTCTTCTTCTTGTTCTTCTTCTCCATCATTATCTTCTTCATCTTCTTC...CGA', 130036: 'AGTCCGAATCTTTTTGCTTGGCTCCTCTGGGGTCACTTCAACCTTATTAATTTC...GGT'}, length=130198)")
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
        self.assertEqual(alignment.query.annotations["splicesites"], (459, 459))
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 263)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({25: 'AGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATG...CAA'}, length=406)")
        self.assertEqual(repr(alignment.target.seq), "Seq({11338: 'TTGTATCTACTTCGTTTTGGCTCTTATAATTTTTTAGTGGCCTAAATTTTCAAA...TGT', 120610: 'TTCTTCTTCTTCTTCTTGTTCTTCTTCTCCATCATTATCTTCTTCATCTTCTTC...CGA', 130036: 'AGTCCGAATCTTTTTGCTTGGCTCCTCTGGGGTCACTTCAACCTTATTAATTTC...GGT'}, length=130198)")
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
        self.assertEqual(alignment.query.annotations["splicesites"], (459, 459))
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_q_multiple_vulgar(self):
        """Test parsing exn_22_q_multiple_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_q_multiple_vulgar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m est2genome comb.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showcigar yes --showvulgar yes")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cdc28p (CDC28) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443482|ref|NC_001134.8|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome II, complete sequence")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 4485)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGAGGGTATTCCAAAGGACCAACCGTTAGGAGCTGATATTGTTAAGAAGTTTATGATGCAACTTTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAAAGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGGAAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAGTGGTACCAAGTCTAGATCCACGCGGTATTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCCAGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATAA")
        self.assertEqual(repr(alignment.target.seq), "Seq({560077: 'ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACA...TAA'}, length=560974)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[560077, 560974], [0, 897]])
            )
        )
        self.assertEqual(alignment[0], "ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGAGGGTATTCCAAAGGACCAACCGTTAGGAGCTGATATTGTTAAGAAGTTTATGATGCAACTTTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAAAGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGGAAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAGTGGTACCAAGTCTAGATCCACGCGGTATTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCCAGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATAA")
        self.assertEqual(alignment[1], "ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGAGGGTATTCCAAAGGACCAACCGTTAGGAGCTGATATTGTTAAGAAGTTTATGATGCAACTTTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAAAGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGGAAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAGTGGTACCAAGTCTAGATCCACGCGGTATTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCCAGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATAA")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cdc28p (CDC28) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443482|ref|NC_001134.8|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome II, complete sequence")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 4485)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGAGGGTATTCCAAAGGACCAACCGTTAGGAGCTGATATTGTTAAGAAGTTTATGATGCAACTTTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAAAGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGGAAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAGTGGTACCAAGTCTAGATCCACGCGGTATTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCCAGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATAA")
        self.assertEqual(repr(alignment.target.seq), "Seq({560077: 'ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACA...TAA'}, length=560974)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[560077, 560974], [0, 897]])
            )
        )
        self.assertEqual(alignment[0], "ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGAGGGTATTCCAAAGGACCAACCGTTAGGAGCTGATATTGTTAAGAAGTTTATGATGCAACTTTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAAAGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGGAAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAGTGGTACCAAGTCTAGATCCACGCGGTATTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCCAGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATAA")
        self.assertEqual(alignment[1], "ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGAGGGTATTCCAAAGGACCAACCGTTAGGAGCTGATATTGTTAAGAAGTTTATGATGCAACTTTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAAAGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGGAAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAGTGGTACCAAGTCTAGATCCACGCGGTATTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCCAGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATAA")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cdc28p (CDC28) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443482|ref|NC_001134.8|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome II, complete sequence")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 4485)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGAGGGTATTCCAAAGGACCAACCGTTAGGAGCTGATATTGTTAAGAAGTTTATGATGCAACTTTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAAAGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGGAAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAGTGGTACCAAGTCTAGATCCACGCGGTATTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCCAGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATAA")
        self.assertEqual(repr(alignment.target.seq), "Seq({560077: 'ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACA...TAA'}, length=560974)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[560077, 560974], [0, 897]])
            )
        )
        self.assertEqual(alignment[0], "ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGAGGGTATTCCAAAGGACCAACCGTTAGGAGCTGATATTGTTAAGAAGTTTATGATGCAACTTTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAAAGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGGAAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAGTGGTACCAAGTCTAGATCCACGCGGTATTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCCAGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATAA")
        self.assertEqual(alignment[1], "ATGAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGAGGGTATTCCAAAGGACCAACCGTTAGGAGCTGATATTGTTAAGAAGTTTATGATGCAACTTTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAAAGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGGAAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAGTGGTACCAAGTCTAGATCCACGCGGTATTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCCAGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATAA")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cdc28p (CDC28) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443753|ref|NC_001148.4|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XVI, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 941)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({2: 'GAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATA...ATA'}, length=896)")
        self.assertEqual(repr(alignment.target.seq), "Seq({492033: 'TAGTACTCTGCAAACCAAGGGTGATGCAGAGCCTGCTTGGCGCTCAGCCTCATA...CTC'}, length=492933)")
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
        self.assertEqual(alignment[0], "GAGC--TCTATGAACAGATTTAAGCAG-TTAGAAAAGCTTGGCAATGGTACGTATGCCACAGTGTACAAGGGACTGAAC--AAAACC---ACAGGGGTATA---T-GTTGCCCTGAAAGAGGTAAAACTGGA---TTCAGAGGAAGGTACACCCTCTACGGCCATCCGTGAGATCTCCCTAATGAAAGAATTGAAACATGAGAACATTGTTAGACTTTATGACGTTATTCACACAGA-G-A-ACAAGTTGACTTTGGTTTTTGAATTCATGGACAACGATTTAAAGAAATACATGGATTCCCGCA--CCGTGGGCAACACACCAAGAGG-GCTAGAACTAAACTTGGTTAAATACTTCCAGTGGCAACTACTGCAAGGG-CTGGCCTTTTGCCATGAAAACAAGATTCTCCACCGTGATTTAAAACCTCAAAACCTATTAATCAACAAGAGA-GGCCAGTTGAAATTGGGTGATTTCGGTCTGGCCCGTGCTTTCGGTATTCCGGTCAACACATTTTCAAGCGAAGTCGTAACGTTGTGGTACCGTGCTCCTGATGTGCTAATGGGTTCTAGGACGTACTCCACATCCATTGATATATGGTCGTGTGGGTGCATTCTTGCGGAAATG-ATAACGGGTAAGCCTTTGTTTCCTGGCACCAACGACGAAGAACAACTGAAATTGATCTTCGACATCATGGGCACTCCTAATGAGTCCCTATGGCCCAGTGTAACAAAGTTACCCAAATACAACCCAAATATCCAGCAACGACCACCAAGAGACCTACGTCAAGTATTGCAACCACACACCAAAGAACCGCTAGACGGGAATCTCATGGATTTCTTACACGGACTCTTGCAACTTAATCCGGATATGAGGCTGAGCGCCAAGCAGGCTCTGCAT-CACCCTTGGTTTGCAGAGTACTA")
        self.assertEqual(alignment[1], "GAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGA---GGGTATTCCAAAGG-ACCA-ACCGTTAGGAGCT-GA--TA---TTG-TTAAGAAGTTTATGATGCAACT-TTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAA-AGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGG-AAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAG---TGGTACCA-AGTC--TAGATCCAC--G-CGGTA------TTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCC-AGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATA")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cdc28p (CDC28) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443753|ref|NC_001148.4|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XVI, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 941)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({2: 'GAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATA...ATA'}, length=896)")
        self.assertEqual(repr(alignment.target.seq), "Seq({492033: 'TAGTACTCTGCAAACCAAGGGTGATGCAGAGCCTGCTTGGCGCTCAGCCTCATA...CTC'}, length=492933)")
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
        self.assertEqual(alignment[0], "GAGC--TCTATGAACAGATTTAAGCAG-TTAGAAAAGCTTGGCAATGGTACGTATGCCACAGTGTACAAGGGACTGAAC--AAAACC---ACAGGGGTATA---T-GTTGCCCTGAAAGAGGTAAAACTGGA---TTCAGAGGAAGGTACACCCTCTACGGCCATCCGTGAGATCTCCCTAATGAAAGAATTGAAACATGAGAACATTGTTAGACTTTATGACGTTATTCACACAGA-G-A-ACAAGTTGACTTTGGTTTTTGAATTCATGGACAACGATTTAAAGAAATACATGGATTCCCGCA--CCGTGGGCAACACACCAAGAGG-GCTAGAACTAAACTTGGTTAAATACTTCCAGTGGCAACTACTGCAAGGG-CTGGCCTTTTGCCATGAAAACAAGATTCTCCACCGTGATTTAAAACCTCAAAACCTATTAATCAACAAGAGA-GGCCAGTTGAAATTGGGTGATTTCGGTCTGGCCCGTGCTTTCGGTATTCCGGTCAACACATTTTCAAGCGAAGTCGTAACGTTGTGGTACCGTGCTCCTGATGTGCTAATGGGTTCTAGGACGTACTCCACATCCATTGATATATGGTCGTGTGGGTGCATTCTTGCGGAAATG-ATAACGGGTAAGCCTTTGTTTCCTGGCACCAACGACGAAGAACAACTGAAATTGATCTTCGACATCATGGGCACTCCTAATGAGTCCCTATGGCCCAGTGTAACAAAGTTACCCAAATACAACCCAAATATCCAGCAACGACCACCAAGAGACCTACGTCAAGTATTGCAACCACACACCAAAGAACCGCTAGACGGGAATCTCATGGATTTCTTACACGGACTCTTGCAACTTAATCCGGATATGAGGCTGAGCGCCAAGCAGGCTCTGCAT-CACCCTTGGTTTGCAGAGTACTA")
        self.assertEqual(alignment[1], "GAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGA---GGGTATTCCAAAGG-ACCA-ACCGTTAGGAGCT-GA--TA---TTG-TTAAGAAGTTTATGATGCAACT-TTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAA-AGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGG-AAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAG---TGGTACCA-AGTC--TAGATCCAC--G-CGGTA------TTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCC-AGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATA")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cdc28p (CDC28) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443753|ref|NC_001148.4|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XVI, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 941)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({2: 'GAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATA...ATA'}, length=896)")
        self.assertEqual(repr(alignment.target.seq), "Seq({492033: 'TAGTACTCTGCAAACCAAGGGTGATGCAGAGCCTGCTTGGCGCTCAGCCTCATA...CTC'}, length=492933)")
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
        self.assertEqual(alignment[0], "GAGC--TCTATGAACAGATTTAAGCAG-TTAGAAAAGCTTGGCAATGGTACGTATGCCACAGTGTACAAGGGACTGAAC--AAAACC---ACAGGGGTATA---T-GTTGCCCTGAAAGAGGTAAAACTGGA---TTCAGAGGAAGGTACACCCTCTACGGCCATCCGTGAGATCTCCCTAATGAAAGAATTGAAACATGAGAACATTGTTAGACTTTATGACGTTATTCACACAGA-G-A-ACAAGTTGACTTTGGTTTTTGAATTCATGGACAACGATTTAAAGAAATACATGGATTCCCGCA--CCGTGGGCAACACACCAAGAGG-GCTAGAACTAAACTTGGTTAAATACTTCCAGTGGCAACTACTGCAAGGG-CTGGCCTTTTGCCATGAAAACAAGATTCTCCACCGTGATTTAAAACCTCAAAACCTATTAATCAACAAGAGA-GGCCAGTTGAAATTGGGTGATTTCGGTCTGGCCCGTGCTTTCGGTATTCCGGTCAACACATTTTCAAGCGAAGTCGTAACGTTGTGGTACCGTGCTCCTGATGTGCTAATGGGTTCTAGGACGTACTCCACATCCATTGATATATGGTCGTGTGGGTGCATTCTTGCGGAAATG-ATAACGGGTAAGCCTTTGTTTCCTGGCACCAACGACGAAGAACAACTGAAATTGATCTTCGACATCATGGGCACTCCTAATGAGTCCCTATGGCCCAGTGTAACAAAGTTACCCAAATACAACCCAAATATCCAGCAACGACCACCAAGAGACCTACGTCAAGTATTGCAACCACACACCAAAGAACCGCTAGACGGGAATCTCATGGATTTCTTACACGGACTCTTGCAACTTAATCCGGATATGAGGCTGAGCGCCAAGCAGGCTCTGCAT-CACCCTTGGTTTGCAGAGTACTA")
        self.assertEqual(alignment[1], "GAGCGGTGAATTAGCAAATTACAAAAGACTTGAGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGACCTGGCCAAGGTCAAAGAGTAGTCGCATTGAAGAAAATAAGACTAGAGAGTGAAGACGAGGGTGTTCCCAGTACAGCCATCAGAGAAATCTCATTATTGAAGGAATTAAAAGACGATAATATTGTCAGATTATACGATATTGTTCACTCTGATGCACACAAGCTATATCTTGTTTTTGAGTTCCTCGATTTGGACCTGAAAAGATATATGGA---GGGTATTCCAAAGG-ACCA-ACCGTTAGGAGCT-GA--TA---TTG-TTAAGAAGTTTATGATGCAACT-TTGTAAGGGTATTGCATACTGCCACTCACACCGTATTCTGCATCGTGATTTAAAACCGCAGAACTTATTGATTAACAA-AGATGGGAATCTAAAACTAGGTGATTTTGGCTTAGCGCGTGCTTTTGGTGTTCCGTTGAGAGCTTACACACATGAAATTGTTACTCTATGGTATAGAGCTCCGGAGGTATTACTGGGTGGAAAACAATATAGTACAGGTGTCGATACATGGTCCATCGGCTGTATATTTGCCGAAATGTGTAACAGG-AAACCAATCTTCAGTGGCGATAGTGAGATCGATCAGATTTTCAAGATATTCAGAGTATTGGGAACGCCGAATGAAGCTATATGGCCAGATATTGTCTACTTGCCTGATTTCAAGCCAAGCTTTCCTCAATGGCGCAGAAAAGACCTATCACAAG---TGGTACCA-AGTC--TAGATCCAC--G-CGGTA------TTGATTTGTTGGACAAACTCCTCGCGTATGACCCTATTAACCGGATTAGCGCC-AGAAGAGCAGCCATCCACCCCTACTTCCAAGAATCATA")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cdc28p (CDC28) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 651)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({34: 'AGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGAC...TCT'}, length=721)")
        self.assertEqual(repr(alignment.target.seq), "Seq({267809: 'AGAAAGTTGGTGAGGGTACTTATGCGGTTGTTTACTTGGGTTGTCAACACTCTA...CCT', 300684: 'AACAAACCGAGATAGAATGGCCAGATATTCTCT'}, length=300717)")
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
        self.assertEqual(alignment.query.annotations["splicesites"], (745,))
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cdc28p (CDC28) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 651)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({34: 'AGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGAC...TCT'}, length=721)")
        self.assertEqual(repr(alignment.target.seq), "Seq({267809: 'AGAAAGTTGGTGAGGGTACTTATGCGGTTGTTTACTTGGGTTGTCAACACTCTA...CCT', 300684: 'AACAAACCGAGATAGAATGGCCAGATATTCTCT'}, length=300717)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
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
                # fmt: off
# flake8: noqa
            )
        )
        self.assertEqual(alignment.query.annotations["splicesites"], (745,))
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296142823|ref|NM_001178508.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cdc28p (CDC28) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 651)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({34: 'AGAAAGTCGGTGAAGGTACATACGGTGTTGTTTATAAAGCGTTAGACTTAAGAC...TCT'}, length=721)")
        self.assertEqual(repr(alignment.target.seq), "Seq({267809: 'AGAAAGTTGGTGAGGGTACTTATGCGGTTGTTTACTTGGGTTGTCAACACTCTA...CCT', 300684: 'AACAAACCGAGATAGAATGGCCAGATATTCTCT'}, length=300717)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
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
                # fmt: off
# flake8: noqa
            )
        )
        self.assertEqual(alignment.query.annotations["splicesites"], (745,))
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 6150)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(repr(alignment.target.seq), "Seq({1318045: 'CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTAT...CAT'}, length=1319275)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        self.assertEqual(alignment[0], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(alignment[1], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 6150)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(repr(alignment.target.seq), "Seq({1318045: 'CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTAT...CAT'}, length=1319275)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        self.assertEqual(alignment[0], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(alignment[1], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 6150)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(repr(alignment.target.seq), "Seq({1318045: 'CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTAT...CAT'}, length=1319275)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates, numpy.array([[1319275, 1318045], [0, 1230]])
            )
        )
        self.assertEqual(alignment[0], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        self.assertEqual(alignment[1], "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAATATTAGCCTTCCTCGATGATCTGCAAGAACAACAGAAAAGGGAAAACGAAAAAGGAACAAGTACAGCAGTTAGCAAGGCTGCAAAGGAATTGCCATCGCCTAATTCAGATGAAAACATGACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGTGAAATGCCAGATTTCAGTCTCGATCTTGTCTCCAATGATAGACAGACTGGTCTAGAAGCTTTAGATTACGACATTCATAACTACTTTCCTCAGCACTCTGAACGCCTGACCGCTGAAAAAATAGATACGTCAGCATGTCAATGTGAAATTGACCAAAAGTATCTTCCATACGAGACAGAAGATGATACTTTATTCCCCAGCGTGCTTCCCCTTGCTGTAGGGAGCCAGTGTAATAATATTTGCAACCGCAAGTGTATCGGGACCAAACCATGTTCAAATAAGGAGATCAAATGCGACTTAATAACAAGCCACCTGTTGAATCAGAAATCTCTAGCTTCGGTGCTTCCGGTGGCTGCTTCTCATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCCATATCGAATGGAAAAGCGTCTTGCTACCACATTCTCGAAGAGATCTCCTCCCTACCAAAATATTCATCGTTGGACATAGATGATTTATGCAGCGAATTAATAATCAAGGCAAAATGTACAGATGACTGCAAAATAGTAGTCAAAGCTCGCGACTTACAGAGTGCTCTGGTTAGACAGCTCCTGTAG")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 439)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAA")
        self.assertEqual(repr(alignment.target.seq), "Seq({85010: 'ATGGTGAACCTCTTCAAGACGGTCAGAATAATCAACAGGATGAAGAAGCAAAAGATGT', 253972: 'TGGCGAGGATAGCGAGCAACCGAAGAAGAAGGGTAGCAAAACTAGCAAAAAGCA...ACA', 350957: 'TGGGGTGATTATATCATTTCTGGATGAGGAATACCTGAAGACCACTCTTCATTA...GTT', 473168: 'AGATGGAAGAATTCTGATAATGCTGTAAAAGAA'}, length=473201)")
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
        self.assertEqual(alignment.query.annotations["splicesites"], (448, 448, 448))
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 439)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAA")
        self.assertEqual(repr(alignment.target.seq), "Seq({85010: 'ATGGTGAACCTCTTCAAGACGGTCAGAATAATCAACAGGATGAAGAAGCAAAAGATGT', 253972: 'TGGCGAGGATAGCGAGCAACCGAAGAAGAAGGGTAGCAAAACTAGCAAAAAGCA...ACA', 350957: 'TGGGGTGATTATATCATTTCTGGATGAGGAATACCTGAAGACCACTCTTCATTA...GTT', 473168: 'AGATGGAAGAATTCTGATAATGCTGTAAAAGAA'}, length=473201)")
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
        self.assertEqual(alignment.query.annotations["splicesites"], (448, 448, 448))
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 439)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ATGGGCAATATCCTTCGGAAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATGGAACACCTAAGAGGAAGGTGGGCAGACCAGGCAGAAAAAGGATTGACTCTGAAGCTAAGAGTAGGAGGACTGCCCAGAATAGGGCAGCTCAACGAGCGTTCCGAGATAGGAAAGAAGCCAAAATGAAGAGTTTGCAAGAGAGGGTAGAGTTACTAGAACAGAAAGATGCGCAGAATAAGACTACCACGGACTTTTTACTATGTTCTTTAAAAAGTTTACTGTCGGAAATTACAAAATATAGAGCTAAGAATTCTGATGATGAAAGAA")
        self.assertEqual(repr(alignment.target.seq), "Seq({85010: 'ATGGTGAACCTCTTCAAGACGGTCAGAATAATCAACAGGATGAAGAAGCAAAAGATGT', 253972: 'TGGCGAGGATAGCGAGCAACCGAAGAAGAAGGGTAGCAAAACTAGCAAAAAGCA...ACA', 350957: 'TGGGGTGATTATATCATTTCTGGATGAGGAATACCTGAAGACCACTCTTCATTA...GTT', 473168: 'AGATGGAAGAATTCTGATAATGCTGTAAAAGAA'}, length=473201)")
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
        self.assertEqual(alignment.query.annotations["splicesites"], (448, 448, 448))
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 263)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({25: 'AGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATG...CAA'}, length=406)")
        self.assertEqual(repr(alignment.target.seq), "Seq({11338: 'TTGTATCTACTTCGTTTTGGCTCTTATAATTTTTTAGTGGCCTAAATTTTCAAA...TGT', 120610: 'TTCTTCTTCTTCTTCTTGTTCTTCTTCTCCATCATTATCTTCTTCATCTTCTTC...CGA', 130036: 'AGTCCGAATCTTTTTGCTTGGCTCCTCTGGGGTCACTTCAACCTTATTAATTTC...GGT'}, length=130198)")
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
        self.assertEqual(alignment.query.annotations["splicesites"], (459, 459))
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 263)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({25: 'AGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATG...CAA'}, length=406)")
        self.assertEqual(repr(alignment.target.seq), "Seq({11338: 'TTGTATCTACTTCGTTTTGGCTCTTATAATTTTTTAGTGGCCTAAATTTTCAAA...TGT', 120610: 'TTCTTCTTCTTCTTCTTGTTCTTCTTCTCCATCATTATCTTCTTCATCTTCTTC...CGA', 130036: 'AGTCCGAATCTTTTTGCTTGGCTCCTCTGGGGTCACTTCAACCTTATTAATTTC...GGT'}, length=130198)")
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
        self.assertEqual(alignment.query.annotations["splicesites"], (459, 459))
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "est2genome")
        self.assertEqual(alignment.score, 263)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({25: 'AGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAATAAAGATG...CAA'}, length=406)")
        self.assertEqual(repr(alignment.target.seq), "Seq({11338: 'TTGTATCTACTTCGTTTTGGCTCTTATAATTTTTTAGTGGCCTAAATTTTCAAA...TGT', 120610: 'TTCTTCTTCTTCTTCTTGTTCTTCTTCTCCATCATTATCTTCTTCATCTTCTTC...CGA', 130036: 'AGTCCGAATCTTTTTGCTTGGCTCCTCTGGGGTCACTTCAACCTTATTAATTTC...GGT'}, length=130198)")
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
        self.assertEqual(alignment.query.annotations["splicesites"], (459, 459))
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_coding2coding_fshifts_cigar(self):
        """Test parsing exn_22_m_cigar_fshifts.exn)."""
        exn_file = os.path.join("Exonerate", "exn_22_m_cigar_fshifts.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m coding2coding c2c_frameshift2.fa scer_cad1.fa --showvulgar no --showalignment yes --bestn 3")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.annotations["model"], "coding2coding")
        self.assertEqual(alignment.score, 213)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAAATTCGTGGGATACTCCCCCCGGAAATCGAACATGCCGTTACCATGACGAAAGTATTAATGGTAGT")
        self.assertEqual(repr(alignment.target.seq), "Seq({465: 'ACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTT...AGT'}, length=630)")
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
        self.assertEqual(alignment[0], "ACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACC-AATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGT")
        self.assertEqual(alignment[1], "ACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAAATTCGTGGGATACTCCCCCCGGAAATCGAACA--TGCCGTTACCAT----GACGAAAGTATTAATGGTAGT")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds:[revcomp]")
        self.assertEqual(alignment.target.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "coding2coding")
        self.assertEqual(alignment.score, 201)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({1: 'CTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTA...GTA'}, length=158)")
        self.assertEqual(repr(alignment.target.seq), "Seq({466: 'CTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTA...GTA'}, length=628)")
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
        self.assertEqual(alignment[0], "TACCATTAATACTTTCGTCACCGATGGTAACGGCACCTGTTCGATTTCCGGGGGGAGAATCCCACGAAT-TGGTTAAACTGGGAGCGTTCCATGAGCCTATGTTCCACATAACTTTCTCATTCTCTTGAGTGTGCGGCTGTACTTCTATACTTGTGTTCACAG")
        self.assertEqual(alignment[1], "TACCATTAATACTTTCGTCATGGT------AACGGCATGTTCGATTTCCGGGGGGAGTATCCCACGAATTTGGTTAAACTGGGAGCGTTCCATGAGCCTATGTTCCACATAACTTTCTCATTCTCTTGAGTGTGCGGCTGTACTTCTATACTTGTGTTCACAG")
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_coding2coding_fshifts_vulgar(self):
        """Test parsing exn_22_o_vulgar_fshifts.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_o_vulgar_fshifts.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m coding2coding c2c_frameshift2.fa scer_cad1.fa --showvulgar no --showalignment yes --bestn 3")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.target.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.annotations["model"], "coding2coding")
        self.assertEqual(alignment.score, 213)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "ACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAAATTCGTGGGATACTCCCCCCGGAAATCGAACATGCCGTTACCATGACGAAAGTATTAATGGTAGT")
        self.assertEqual(repr(alignment.target.seq), "Seq({465: 'ACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTT...AGT'}, length=630)")
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
        self.assertEqual(alignment[0], "ACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACC-AATTCGTGGGATTCTCCCCCCGGAAATCGAACAGGTGCCGTTACCATCGGTGACGAAAGTATTAATGGTAGT")
        self.assertEqual(alignment[1], "ACTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTATGTGGAACATAGGCTCATGGAACGCTCCCAGTTTAACCAAATTCGTGGGATACTCCCCCCGGAAATCGAACA--TGCCGTTACCAT----GACGAAAGTATTAATGGTAGT")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.query.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds:[revcomp]")
        self.assertEqual(alignment.target.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "coding2coding")
        self.assertEqual(alignment.score, 201)
        self.assertGreater(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({1: 'CTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTA...GTA'}, length=158)")
        self.assertEqual(repr(alignment.target.seq), "Seq({466: 'CTGTGAACACAAGTATAGAAGTACAGCCGCACACTCAAGAGAATGAGAAAGTTA...GTA'}, length=628)")
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
        self.assertEqual(alignment[0], "TACCATTAATACTTTCGTCACCGATGGTAACGGCACCTGTTCGATTTCCGGGGGGAGAATCCCACGAAT-TGGTTAAACTGGGAGCGTTCCATGAGCCTATGTTCCACATAACTTTCTCATTCTCTTGAGTGTGCGGCTGTACTTCTATACTTGTGTTCACAG")
        self.assertEqual(alignment[1], "TACCATTAATACTTTCGTCATGGT------AACGGCATGTTCGATTTCCGGGGGGAGTATCCCACGAATTTGGTTAAACTGGGAGCGTTCCATGAGCCTATGTTCCACATAACTTTCTCATTCTCTTGAGTGTGCGGCTGTACTTCTATACTTGTGTTCACAG")
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_protein2dna_cigar(self):
        """Test parsing exonerate output (exn_22_m_protein2dna_cigar.exn)."""
        exn_file = os.path.join("Exonerate", "exn_22_m_protein2dna_cigar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m protein2dna ../scer_cad1_prot.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showcigar no --showvulgar no")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.query.description, "AP-1-like transcription activator YAP2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) GN=CAD1 PE=1 SV=2")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence")
        self.assertEqual(alignment.annotations["model"], "protein2dna:local")
        self.assertEqual(alignment.score, 2105)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[0, 409], [1319275, 1318048]])
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.query.description, "AP-1-like transcription activator YAP2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) GN=CAD1 PE=1 SV=2")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence")
        self.assertEqual(alignment.annotations["model"], "protein2dna:local")
        self.assertEqual(alignment.score, 205)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal( alignment.coordinates, numpy.array([[0, 0], [0, 0]])
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.query.description, "AP-1-like transcription activator YAP2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) GN=CAD1 PE=1 SV=2")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence")
        self.assertEqual(alignment.annotations["model"], "protein2dna:local")
        self.assertEqual(alignment.score, 116)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal( alignment.coordinates, numpy.array([[0, 0], [0, 0]])
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_protein2dna_vulgar(self):
        """Test parsing exonerate output (exn_22_m_protein2dna_vulgar.exn)."""
        exn_file = os.path.join("Exonerate", "exn_22_m_protein2dna_vulgar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m protein2dna ../scer_cad1_prot.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showcigar no --showvulgar no")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.query.description, "AP-1-like transcription activator YAP2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) GN=CAD1 PE=1 SV=2")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence")
        self.assertEqual(alignment.annotations["model"], "protein2dna:local")
        self.assertEqual(alignment.score, 2105)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[0, 409], [1319275, 1318048]])
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.query.description, "AP-1-like transcription activator YAP2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) GN=CAD1 PE=1 SV=2")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence")
        self.assertEqual(alignment.annotations["model"], "protein2dna:local")
        self.assertEqual(alignment.score, 205)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal( alignment.coordinates, numpy.array([[0, 0], [0, 0]])
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.query.description, "AP-1-like transcription activator YAP2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) GN=CAD1 PE=1 SV=2")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence")
        self.assertEqual(alignment.annotations["model"], "protein2dna:local")
        self.assertEqual(alignment.score, 116)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertTrue(
            numpy.array_equal( alignment.coordinates, numpy.array([[0, 0], [0, 0]])
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_protein2dna_fshifts_cigar(self):
        """Test parsing exonerate output (exn_22_o_cigar_fshifts2.exn)."""
        exn_file = os.path.join("Exonerate", "exn_22_o_cigar_fshifts2.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m protein2dna scer_cad1_prot.fa scer_cad1frameshift.fa --showvulgar no --bestn 3")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.query.description, "AP-1-like transcription activator YAP2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) GN=CAD1 PE=1 SV=2")
        self.assertEqual(alignment.target.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.annotations["model"], "protein2dna:local")
        self.assertEqual(alignment.score, 367)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({330: 'HTKTIRTQSEAIEHISSAISNGKASCYHILEEISSLPKYSSLDIDDLCSELIIKAKCTDDCKIVVKARDLQSALVRQLL'})", length=409)
        self.assertEqual(repr(alignment.target.seq), "Seq({216: 'CATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCC...CTG'}, length=455)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[330, 373, 373, 409], [216, 345, 347, 455]])
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.query.description, "AP-1-like transcription activator YAP2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) GN=CAD1 PE=1 SV=2")
        self.assertEqual(alignment.target.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.annotations["model"], "protein2dna:local")
        self.assertEqual(alignment.score, 322)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({6: 'KGQQIYLAGDMKKQMLLNKDGTPKRKVGRPGRKRIDSEAKSRRTAQNRAAQRAF...MKS'}, length=70)")
        self.assertEqual(repr(alignment.target.seq), "Seq({16: 'AAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAAT...AGT'}, length=208)")
        self.assertTrue(
            numpy.array_equal( alignment.coordinates, numpy.array([[6, 70], [16, 208]])
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_protein2dna_fshifts_vulgar(self):
        """Test parsing exonerate output (exn_22_o_vulgar_fshifts2.exn)."""
        exn_file = os.path.join("Exonerate", "exn_22_o_vulgar_fshifts2.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m protein2dna scer_cad1_prot.fa scer_cad1frameshift.fa --showvulgar no --bestn 3")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.query.description, "AP-1-like transcription activator YAP2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) GN=CAD1 PE=1 SV=2")
        self.assertEqual(alignment.target.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.annotations["model"], "protein2dna:local")
        self.assertEqual(alignment.score, 367)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({330: 'HTKTIRTQSEAIEHISSAISNGKASCYHILEEISSLPKYSSLDIDDLCSELIIKAKCTDDCKIVVKARDLQSALVRQLL'})", length=409)
        self.assertEqual(repr(alignment.target.seq), "Seq({216: 'CATACTAAAACAATTCGAACCCAATCTGAAGCAATTGAACACATTAGCAGCGCC...CTG'}, length=455)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[330, 373, 373, 409], [216, 345, 347, 455]])
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.query.description, "AP-1-like transcription activator YAP2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) GN=CAD1 PE=1 SV=2")
        self.assertEqual(alignment.target.id, "gi|296143771|ref|NM_001180731.1|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c Cad1p (CAD1) mRNA, complete cds")
        self.assertEqual(alignment.annotations["model"], "protein2dna:local")
        self.assertEqual(alignment.score, 322)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({6: 'KGQQIYLAGDMKKQMLLNKDGTPKRKVGRPGRKRIDSEAKSRRTAQNRAAQRAF...MKS'}, length=70)")
        self.assertEqual(repr(alignment.target.seq), "Seq({16: 'AAAGGTCAGCAAATATATTTAGCAGGTGACATGAAGAAGCAAATGTTGCTAAAT...AGT'}, length=208)")
        self.assertTrue(
            numpy.array_equal( alignment.coordinates, numpy.array([[6, 70], [16, 208]])
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_protein2genome_cigar(self):
        """Test parsing exn_22_m_protein2genome_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_protein2genome_cigar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m protein2genome ../scer_cad1_prot.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showcigar no --showvulgar no")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.query.description, "AP-1-like transcription activator YAP2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) GN=CAD1 PE=1 SV=2")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "protein2genome:local")
        self.assertEqual(alignment.score, 2105)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "Seq('MGNILRKGQQIYLAGDMKKQMLLNKDGTPKRKVGRPGRKRIDSEAKSRRTAQNRAAQRAFRDRKEAKMKSLQERVELLEQKDAQNKTTTDFLLCSLKSLLSEITKYRAKNSDDERILAFLDDLQEQQKRENEKGTSTAVSKAAKELPSPNSDENMTVNTSIEVQPHTQENEKVMWNIGSWNAPSLTNSWDSPPGNRTGAVTIGDESINGSEMPDFSLDLVSNDRQTGLEALDYDIHNYFPQHSERLTAEKIDTSACQCEIDQKYLPYETEDDTLFPSVLPLAVGSQCNNICNRKCIGTKPCSNKEIKCDLITSHLLNQKSLASVLPVAASHTKTIRTQSEAIEHISSAISNGKASCYHILEEISSLPKYSSLDIDDLCSELIIKAKCTDDCKIVVKARDLQSALVRQLL')", length=409)
        self.assertEqual(repr(alignment.target.seq), "Seq({1318045: 'CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTAT...CAT'}, length=1319275)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[1319275, 1318048], [0, 409]])
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.query.description, "AP-1-like transcription activator YAP2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) GN=CAD1 PE=1 SV=2")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence")
        self.assertEqual(alignment.annotations["model"], "protein2genome:local")
        self.assertEqual(alignment.score, 205)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({28: 'PKRKVGRPGRKRIDSEAKSRRTAQNRAAQRAFRDRKEAKMKSLQERVELLEQKD...AFL'}, length=120)")
        self.assertEqual(repr(alignment.target.seq), "Seq({253991: 'CCGAAGAAGAAGGGTAGCAAAACTAGCAAAAAGCAAGATTTGGATCCTGAAACT...TTA'}, length=254270)")
        self.assertTrue(
            numpy.array_equal( alignment.coordinates, numpy.array([[28, 120],[253991, 254270]])
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.query.description, "AP-1-like transcription activator YAP2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) GN=CAD1 PE=1 SV=2")
        self.assertEqual(alignment.target.id, "gi|330443590|ref|NC_001140.6|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome VIII, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "protein2genome:local")
        self.assertEqual(alignment.score, 122)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({37: 'RKRIDSEAKSRRTAQNRAAQRAFRDRKEAKMKSLQERVELLEQKDAQNKTTTDF...LQE'}, length=125)")
        self.assertEqual(repr(alignment.target.seq), "Seq({68450: 'ATAAACAAGATACTGAACAGGGACCCCCAGTTCATGTCGAATTCCAGCTTTCAT...GAA', 84531: 'AATGAGAATGTTCCAGATGACTCTAAGGCAAAGAAAAAGGCTCAAAACAGAGCC...TGT'}, length=84646)")
        self.assertTrue(
            numpy.array_equal( alignment.coordinates,
            numpy.array([[   37,    75,    75,    75,   125],
                         [84646, 84533, 68601, 68600, 68450]])
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_m_protein2genome_vulgar(self):
        """Test parsing exn_22_m_protein2genome_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_22_m_protein2genome_vulgar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m protein2genome ../scer_cad1_prot.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showcigar no --showvulgar no")
        self.assertEqual(alignments.hostname, "blackbriar")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.query.description, "AP-1-like transcription activator YAP2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) GN=CAD1 PE=1 SV=2")
        self.assertEqual(alignment.target.id, "gi|330443520|ref|NC_001136.10|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome IV, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "protein2genome:local")
        self.assertEqual(alignment.score, 2105)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(alignment.query.seq, "Seq('MGNILRKGQQIYLAGDMKKQMLLNKDGTPKRKVGRPGRKRIDSEAKSRRTAQNRAAQRAFRDRKEAKMKSLQERVELLEQKDAQNKTTTDFLLCSLKSLLSEITKYRAKNSDDERILAFLDDLQEQQKRENEKGTSTAVSKAAKELPSPNSDENMTVNTSIEVQPHTQENEKVMWNIGSWNAPSLTNSWDSPPGNRTGAVTIGDESINGSEMPDFSLDLVSNDRQTGLEALDYDIHNYFPQHSERLTAEKIDTSACQCEIDQKYLPYETEDDTLFPSVLPLAVGSQCNNICNRKCIGTKPCSNKEIKCDLITSHLLNQKSLASVLPVAASHTKTIRTQSEAIEHISSAISNGKASCYHILEEISSLPKYSSLDIDDLCSELIIKAKCTDDCKIVVKARDLQSALVRQLL')", length=409)
        self.assertEqual(repr(alignment.target.seq), "Seq({1318045: 'CTACAGGAGCTGTCTAACCAGAGCACTCTGTAAGTCGCGAGCTTTGACTACTAT...CAT'}, length=1319275)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[1319275, 1318048], [0, 409]])
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.query.description, "AP-1-like transcription activator YAP2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) GN=CAD1 PE=1 SV=2")
        self.assertEqual(alignment.target.id, "gi|330443688|ref|NC_001145.3|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome XIII, complete sequence")
        self.assertEqual(alignment.annotations["model"], "protein2genome:local")
        self.assertEqual(alignment.score, 205)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertLess(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({28: 'PKRKVGRPGRKRIDSEAKSRRTAQNRAAQRAFRDRKEAKMKSLQERVELLEQKD...AFL'}, length=120)")
        self.assertEqual(repr(alignment.target.seq), "Seq({253991: 'CCGAAGAAGAAGGGTAGCAAAACTAGCAAAAAGCAAGATTTGGATCCTGAAACT...TTA'}, length=254270)")
        self.assertTrue(
            numpy.array_equal( alignment.coordinates, numpy.array([[28, 120],[253991, 254270]])
            )
        )
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "sp|P24813|YAP2_YEAST")
        self.assertEqual(alignment.query.description, "AP-1-like transcription activator YAP2 OS=Saccharomyces cerevisiae (strain ATCC 204508 / S288c) GN=CAD1 PE=1 SV=2")
        self.assertEqual(alignment.target.id, "gi|330443590|ref|NC_001140.6|")
        self.assertEqual(alignment.target.description, "Saccharomyces cerevisiae S288c chromosome VIII, complete sequence:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "protein2genome:local")
        self.assertEqual(alignment.score, 122)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({37: 'RKRIDSEAKSRRTAQNRAAQRAFRDRKEAKMKSLQERVELLEQKDAQNKTTTDF...LQE'}, length=125)")
        self.assertEqual(repr(alignment.target.seq), "Seq({68450: 'ATAAACAAGATACTGAACAGGGACCCCCAGTTCATGTCGAATTCCAGCTTTCAT...GAA', 84531: 'AATGAGAATGTTCCAGATGACTCTAAGGCAAAGAAAAAGGCTCAAAACAGAGCC...TGT'}, length=84646)")
        self.assertTrue(
            numpy.array_equal( alignment.coordinates,
            numpy.array([[   37,    75,    75,    75,   125],
                         [84646, 84533, 68601, 68600, 68450]])
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_24_m_protein2genome_revcomp_fshifts_cigar(self):
        """Test parsing exn_24_m_protein2genome_revcomp_fshifts_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_24_m_protein2genome_revcomp_fshifts_cigar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m protein2genome gene026_baits.fasta gene026_contigs.fasta --showcigar no --showvulgar no --bestn 2 --refine full")
        self.assertEqual(alignments.hostname, "RBGs-MacBook-Air.local")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "Morus-gene026")
        self.assertEqual(alignment.target.id, "NODE_2_length_1708_cov_48.590765")
        self.assertEqual(alignment.target.description, "SPAdes NODE_2 contig:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "protein2genome:local")
        self.assertEqual(alignment.score, 1308)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({69: 'PESPWTCSPLQTPSPSLLYHCIASLHRHDGTIHSIAVSNGVVFTGSDSRRIRAW...DEE'}, length=441)")
        self.assertEqual(repr(alignment.target.seq), "Seq({331: 'TTCTTGATCAGGCAAGATTTTGACTCTCCAGACCTTCAAAGTCTGGTCCAAGCT...TGG'}, length=1416)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                array([[1416, 1380, 1374, 1125, 1125, 1047, 1047,  744,  744,
                         450,  450, 448,  331],
                       [  69,   81,   81,  164,  169,  195,  196,  297,  300,
                         398,  402, 402,  441]])
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_24_m_protein2genome_revcomp_fshifts_vulgar(self):
        """Test parsing exn_24_m_protein2genome_revcomp_fshifts_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_24_m_protein2genome_revcomp_fshifts_vulgar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m protein2genome gene026_baits.fasta gene026_contigs.fasta --showcigar no --showvulgar no --bestn 2 --refine full")
        self.assertEqual(alignments.hostname, "RBGs-MacBook-Air.local")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "Morus-gene026")
        self.assertEqual(alignment.target.id, "NODE_2_length_1708_cov_48.590765")
        self.assertEqual(alignment.target.description, "SPAdes NODE_2 contig:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "protein2genome:local")
        self.assertEqual(alignment.score, 1308)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({69: 'PESPWTCSPLQTPSPSLLYHCIASLHRHDGTIHSIAVSNGVVFTGSDSRRIRAW...DEE'}, length=441)")
        self.assertEqual(repr(alignment.target.seq), "Seq({331: 'TTCTTGATCAGGCAAGATTTTGACTCTCCAGACCTTCAAAGTCTGGTCCAAGCT...TGG'}, length=1416)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                array([[1416, 1380, 1374, 1125, 1125, 1047, 1047,  744,  744,
                         450,  450, 448,  331],
                       [  69,   81,   81,  164,  169,  195,  196,  297,  300,
                         398,  402, 402,  441]])
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_24_protein2genome_met_intron_cigar(self):
        """Test parsing exn_24_m_protein2genome_met_intron_cigar.exn."""
        exn_file = os.path.join("Exonerate", "exn_24_m_protein2genome_met_intron_cigar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m protein2genome gene001_baits.fasta gene001_contigs.fasta --showcigar no --showvulgar no --bestn 1 --refine full")
        self.assertEqual(alignments.hostname, "RBGs-MacBook-Air.local")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "Morus-gene001")
        self.assertEqual(alignment.target.id, "NODE_1_length_2817_cov_100.387732")
        self.assertEqual(alignment.target.description, "SPAdes contig NODE_1:[revcomp]")
        self.assertEqual(alignment.annotations["model"], "protein2genome:local")
        self.assertEqual(alignment.score, 1958)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({48: 'MVQTPLHVSAGNNRADIVKFLLEFPGPEKVELEAKNMYGETPLHMAAKNGCNDA...SQ*'}, length=482)")
        self.assertEqual(repr(alignment.target.seq), "Seq({388: 'TCATTGCGTCAAAAGCTGAATTCCTGCTTCTAAATCCTCTAGGGTAATTGTACA...TCT', 1056: 'ACCTTCCTCCTAGTTTTTGGTCCAGTATGACCAACAAATTCCCCAACCAAATCA...CCT', 1418: 'ACCTGTTCCAGGGTTGCCAAGGAAGGCCATATGAGGCGGTCTTCGGGTACCAAC...CCT', 1808: 'ACATTGTCCTCTGCACTGCAATCAGCATTGTACTCGAGCAATGTCTTTACAGTT...TCT', 2028: 'ACATTTGCTTTGGCTTCAATAAAAGCACCATGAGCAAGAAGCACCCGAGCAGCA...ACT', 2279: 'ACCATGTTCTTGGCTTCCAACTCAACCTTCTCTGGCCCTGGAAACTCAAGAAGA...CAT'}, length=2392)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                array([[2392, 2281, 2129, 2030, 1921, 1810, 1724, 1420, 1198,
                        1196, 1058,  925,  388],
                       [  48,   85,   85,  118,  118,  155,  155,  257,  257,
                         257,  303,  303,  482]])
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_24_protein2genome_met_intron_vulgar(self):
        """Test parsing exn_24_m_protein2genome_met_intron_vulgar.exn."""
        exn_file = os.path.join("Exonerate", "exn_24_m_protein2genome_met_intron_vulgar.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m protein2genome gene001_baits.fasta gene001_contigs.fasta --showalignment no --showcigar no --showvulgar yes --bestn 1 --refine full")
        self.assertEqual(alignments.hostname, "Michiels-MacBook-Pro.local")
        alignment = next(alignments)
        self.assertEqual(alignment.query.id, "Morus-gene001")
        self.assertEqual(alignment.target.id, "NODE_1_length_2817_cov_100.387732")
        self.assertEqual(alignment.score, 1978)
        self.assertLess(alignment.coordinates[0, 0], alignment.coordinates[0, -1])
        self.assertGreater(alignment.coordinates[1, 0], alignment.coordinates[1, -1])
        self.assertEqual(repr(alignment.query.seq), "Seq({48: 'MVQTPLHVSAGNNRADIVKFLLEFPGPEKVELEAKNMYGETPLHMAAKNGCNDA...SQ*'}, length=482)")
        self.assertEqual(repr(alignment.target.seq), "Seq({388: 'TCATTGCGTCAAAAGCTGAATTCCTGCTTCTAAATCCTCTAGGGTAATTGTACA...TCT', 1056: 'ACCTTCCTCCTAGTTTTTGGTCCAGTATGACCAACAAATTCCCCAACCAAATCA...CCT', 1418: 'ACCTGTTCCAGGGTTGCCAAGGAAGGCCATATGAGGCGGTCTTCGGGTACCAAC...CCT', 1808: 'ACATTGTCCTCTGCACTGCAATCAGCATTGTACTCGAGCAATGTCTTTACAGTT...TCT', 2028: 'ACATTTGCTTTGGCTTCAATAAAAGCACCATGAGCAAGAAGCACCCGAGCAGCA...ACT', 2279: 'ACCATGTTCTTGGCTTCCAACTCAACCTTCTCTGGCCCTGGAAACTCAAGAAGA...CAT'}, length=2392)")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                array([[2392, 2281, 2129, 2030, 1921, 1810, 1724, 1420, 1198,
                        1196, 1058,  925,  388],
                       [  48,   85,   85,  118,  118,  155,  155,  257,  257,
                         257,  303,  303,  482]])
            )
        )
        with self.assertRaises(StopIteration):
            next(alignments)

    def test_exn_22_q_none(self):
        """Test parsing exonerate output (exn_22_q_none.exn)."""
        exn_file = os.path.join("Exonerate", "exn_22_q_none.exn")
        alignments = exonerate.AlignmentIterator(exn_file)
        self.assertEqual(alignments.program, "exonerate")
        self.assertEqual(alignments.commandline, "exonerate -m est2genome none.fa /media/Waterloo/Downloads/genomes/scer_s288c/scer_s288c.fa --bestn 3 --showcigar yes --showvulgar yes")
        self.assertEqual(alignments.hostname, "blackbriar")
        self.assertRaises(StopIteration, next, alignments)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
