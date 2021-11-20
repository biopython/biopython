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


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
