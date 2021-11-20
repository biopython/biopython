# Copyright 2021 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Tests for Align.maf module."""
import unittest


from Bio.Align import maf


try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install numpy if you want to use Bio.Align.maf."
    ) from None


class TestAlignIO_reading(unittest.TestCase):
    def test_reading_bundle_without_target(self):
        """Test parsing bundle_without_target.maf."""
        path = "MAF/bundle_without_target.maf"
        alignments = maf.AlignmentIterator(path)
        self.assertEqual(alignments.metadata["version"], "1")
        self.assertEqual(alignments.metadata["scoring"], "autoMZ.v1")
        alignment = next(alignments)
        self.assertRaises(StopIteration, next, alignments)
        self.assertEqual(alignment.score, 6441)
        self.assertEqual(len(alignment.sequences), 2)
        self.assertEqual(alignment.sequences[0].id, "mm8.chr10")
        self.assertEqual(alignment.sequences[1].id, "oryCun1.scaffold_133159")
        self.assertEqual(len(alignment.sequences[0].seq), 129993255)
        self.assertEqual(len(alignment.sequences[1].seq), 13221)
        self.assertEqual(alignment.sequences[0].seq[3009319: 3009319+162], "TCATAGGTATTTATTTTTAAATATGGTTTGCTTTATGGCTAGAACACACCGATTACTTAAAATAGGATTAACCCCCATACACTTTAAAAATGATTAAACAACATTTCTGCTGCTCGCTCACATTCTTCATAGAAGATGACATAATGTATTTTCCTTTTGGTT")
        self.assertEqual(alignment.sequences[1].seq[11087: 11087+164], "TCACAGATATTTACTATTAAATATGGTTTGTTATATGGTTACGGTTCATAGGTTACTTGGAATTGGATTAACCTTCTTATTCATTGCAGAATTGGTTACACTGTGTTCTTGACCTTTGCTTGTTTTCTCCATGGAAACTGATGTCAAATACTTTCCCTTTGGTT")
        self.assertEqual(alignment[0], "TCATAGGTATTTATTTTTAAATATGGTTTGCTTTATGGCTAGAACACACCGATTACTTAAAATAGGATTAACC--CCCATACACTTTAAAAATGATTAAACAACATTTCTGCTGCTCGCTCACATTCTTCATAGAAGATGACATAATGTATTTTCCTTTTGGTT")
        self.assertEqual(alignment[1], "TCACAGATATTTACTATTAAATATGGTTTGTTATATGGTTACGGTTCATAGGTTACTTGGAATTGGATTAACCTTCTTATTCATTGCAGAATTGGTTACACTGTGTTCTTGACCTTTGCTTGTTTTCTCCATGGAAACTGATGTCAAATACTTTCCCTTTGGTT")
        self.assertEqual(alignment.column_annotations["oryCun1.scaffold_133159"], "99569899999998999999999999999999999999999999999999999999999999999999999757878999975999999999999999979999999999997899999999999997997999999869999996999988997997999999")
        self.assertEqual(alignment.sequences[1].annotations['leftStatus'], 'N')
        self.assertEqual(alignment.sequences[1].annotations['leftCount'], 0)
        self.assertEqual(alignment.sequences[1].annotations['rightStatus'], 'N')
        self.assertEqual(alignment.sequences[1].annotations['rightCount'], 0)
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([[3009319, 3009392, 3009392, 3009481],
                             [  11087,   11160,   11162,   11251],
                            ])
                )
            )

    def test_reading_bug2453(self):
        """Test parsing bug2453.maf."""
        path = "MAF/bug2453.maf"
        alignments = maf.AlignmentIterator(path)
        self.assertEqual(alignments.metadata["version"], "1")
        self.assertEqual(alignments.metadata["scoring"], "tba.v8")
        self.assertEqual(alignments.metadata["comments"], [
                 "tba.v8 (((human chimp) baboon) (mouse rat))",
                 "multiz.v7",
                 "maf_project.v5 _tba_right.maf3 mouse _tba_C",
                 "single_cov2.v4 single_cov2 /dev/stdin",
         ])
        alignment = next(alignments)
        self.assertEqual(alignment.score, 23262)
        self.assertEqual(len(alignment.sequences), 5)
        self.assertEqual(alignment.sequences[0].id, "hg16.chr7")
        self.assertEqual(len(alignment.sequences[0]), 158545518)
        self.assertEqual(alignment.sequences[0].seq[27578828:27578828+38], "AAAGGGAATGTTAACCAAATGAATTGTCTCTTACGGTG")
        self.assertEqual(alignment[0], "AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG")
        self.assertEqual(alignment.sequences[1].id, "panTro1.chr6")
        self.assertEqual(len(alignment.sequences[1]), 161576975)
        self.assertEqual(alignment.sequences[1].seq[28741140:28741140+38], "AAAGGGAATGTTAACCAAATGAATTGTCTCTTACGGTG")
        self.assertEqual(alignment[1], "AAA-GGGAATGTTAACCAAATGA---ATTGTCTCTTACGGTG")
        self.assertEqual(alignment.sequences[2].id, "baboon")
        self.assertEqual(len(alignment.sequences[2]), 4622798)
        self.assertEqual(alignment.sequences[2].seq[116834:116834+38], "AAAGGGAATGTTAACCAAATGAGTTGTCTCTTATGGTG")
        self.assertEqual(alignment[2], "AAA-GGGAATGTTAACCAAATGA---GTTGTCTCTTATGGTG")
        self.assertEqual(alignment.sequences[3].id, "mm4.chr6")
        self.assertEqual(len(alignment.sequences[3]), 151104725)
        self.assertEqual(alignment.sequences[3].seq[53215344:53215344+38], "AATGGGAATGTTAAGCAAACGAATTGTCTCTCAGTGTG")
        self.assertEqual(alignment[3], "-AATGGGAATGTTAAGCAAACGA---ATTGTCTCTCAGTGTG")
        self.assertEqual(alignment.sequences[4].id, "rn3.chr4")
        self.assertEqual(len(alignment.sequences[4]), 187371129)
        self.assertEqual(alignment.sequences[4].seq[81344243:81344243+40], "AAGGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG")
        self.assertEqual(alignment[4], "-AA-GGGGATGCTAAGCCAATGAGTTGTTGTCTCTCAATGTG")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([
        [27578828, 27578829, 27578831, 27578831, 27578850, 27578850, 27578866],
        [28741140, 28741141, 28741143, 28741143, 28741162, 28741162, 28741178],
        [  116834,   116835,   116837,   116837,   116856,   116856,   116872],
        [53215344, 53215344, 53215346, 53215347, 53215366, 53215366, 53215382],
        [81344243, 81344243, 81344245, 81344245, 81344264, 81344267, 81344283],
                            ])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 5062.0 )
        self.assertEqual(len(alignment.sequences), 5)
        self.assertEqual(alignment.sequences[0].id, "hg16.chr7")
        self.assertEqual(len(alignment.sequences[0]), 158545518)
        self.assertEqual(alignment.sequences[0].seq[27699739:27699739+6], "TAAAGA")
        self.assertEqual(alignment[0], "TAAAGA")
        self.assertEqual(alignment.sequences[1].id, "panTro1.chr6")
        self.assertEqual(len(alignment.sequences[1]), 161576975)
        self.assertEqual(alignment.sequences[1].seq[28862317:28862317+6], "TAAAGA")
        self.assertEqual(alignment[1], "TAAAGA")
        self.assertEqual(alignment.sequences[2].id, "baboon")
        self.assertEqual(len(alignment.sequences[2]), 4622798)
        self.assertEqual(alignment.sequences[2].seq[241163:241163+6], "TAAAGA")
        self.assertEqual(alignment[2], "TAAAGA")
        self.assertEqual(alignment.sequences[3].id, "mm4.chr6")
        self.assertEqual(len(alignment.sequences[3]), 151104725)
        self.assertEqual(alignment.sequences[3].seq[53303881:53303881+6], "TAAAGA")
        self.assertEqual(alignment[3], "TAAAGA")
        self.assertEqual(alignment.sequences[4].id, "rn3.chr4")
        self.assertEqual(len(alignment.sequences[4]), 187371129)
        self.assertEqual(alignment.sequences[4].seq[81444246:81444246+6], "taagga")
        self.assertEqual(alignment[4], "taagga")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([
                             [27699739, 27699745],
                             [28862317, 28862323],
                             [  241163,   241169],
                             [53303881, 53303887],
                             [81444246, 81444252],
                            ])
                )
            )
        alignment = next(alignments)
        self.assertEqual(alignment.score, 6636.0)
        self.assertEqual(len(alignment.sequences), 4)
        self.assertEqual(alignment.sequences[0].id, "hg16.chr7")
        self.assertEqual(len(alignment.sequences[0]), 158545518)
        self.assertEqual(alignment.sequences[0].seq[27707221:27707221+13], "gcagctgaaaaca")
        self.assertEqual(alignment[0], "gcagctgaaaaca")
        self.assertEqual(alignment.sequences[1].id, "panTro1.chr6")
        self.assertEqual(len(alignment.sequences[1]), 161576975)
        self.assertEqual(alignment.sequences[1].seq[28869787:28869787+13], "gcagctgaaaaca")
        self.assertEqual(alignment[1], "gcagctgaaaaca")
        self.assertEqual(alignment.sequences[2].id, "baboon")
        self.assertEqual(len(alignment.sequences[2]), 4622798)
        self.assertEqual(alignment.sequences[2].seq[249182:249182+13], "gcagctgaaaaca")
        self.assertEqual(alignment[2], "gcagctgaaaaca")
        self.assertEqual(alignment.sequences[3].id, "mm4.chr6")
        self.assertEqual(len(alignment.sequences[3]), 151104725)
        self.assertEqual(alignment.sequences[3].seq[53310102:53310102+13], "ACAGCTGAAAATA")
        self.assertEqual(alignment[3], "ACAGCTGAAAATA")
        self.assertTrue(
            numpy.array_equal(
                alignment.coordinates,
                numpy.array([
                    [27707221, 27707234],
                    [28869787, 28869800],
                    [  249182,   249195],
                    [53310102, 53310115],
                            ])
                )
            )
        self.assertRaises(StopIteration, next, alignments)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
