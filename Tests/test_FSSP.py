# Copyright 2001 by Iddo Friedberg.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Tests for the FSSP module."""

import os
import unittest

from Bio import FSSP
from Bio.FSSP import FSSPTools


class TestGeo(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        path = os.path.join('FSSP', '1cnv.fssp')
        handle = open(path)
        cls.head_rec, cls.sum_rec, cls.align_rec = FSSP.read_fssp(handle)
        handle.close()

    def test_attributes(self):
        self.assertEqual(self.head_rec.author, ['M.Hennig'])
        self.assertEqual(self.head_rec.compnd, ['concanavalin', 'b'])
        self.assertEqual(self.head_rec.database, 2645)
        self.assertEqual(self.head_rec.header, "SEED PROTEIN")
        self.assertEqual(self.head_rec.nalign, 214)
        self.assertEqual(self.head_rec.pdbid, '1cnv')
        self.assertEqual(self.head_rec.seqlength, 283)
        self.assertEqual(self.head_rec.source, '(canavalia ensiformis) jack bean')

    def test_alignment(self):
        self.assertEqual(len(self.sum_rec), self.head_rec.nalign)
        alignment = FSSPTools.mult_align(self.sum_rec, self.align_rec)
        name_list = ['2hvm0', '1hvq0', '1nar0', '2ebn0']
        sum_newnames, align_newnames = FSSPTools.name_filter(self.sum_rec, self.align_rec,
                                                             name_list)
        self.assertEqual(len(sum_newnames), 4)
        line = """\
   2: 1cnv   2hvm   39.2  1.7  270   273   42      0      0    10 S    hevamine (chitinaseLYSOZYME) 
"""  # noqa : W291
        self.assertEqual(str(sum_newnames[2]), line)
        line = """\
   3: 1cnv   1hvq   39.0  1.7  271   273   41      0      0    10 S    hevamine a 
"""  # noqa : W291
        self.assertEqual(str(sum_newnames[3]), line)
        line = """\
   5: 1cnv   1nar   20.0  3.1  246   289   13      0      0    27 S    Narbonin 
"""  # noqa : W291
        self.assertEqual(str(sum_newnames[5]), line)
        line = """\
  11: 1cnv   2ebn   16.9  3.0  215   285   13      0      0    25 S    Endo-beta-n-acetylglucosaminidase f1 (endoglycosidase f
"""  # noqa : W291
        self.assertEqual(str(sum_newnames[11]), line)
        new_dict = align_newnames['0P168'].pos_align_dict
        self.assertEqual(len(new_dict), 4)
        self.assertEqual(str(new_dict[2]), 'Ps')
        self.assertEqual(str(new_dict[3]), 'Ps')
        self.assertEqual(str(new_dict[5]), '..')
        self.assertEqual(str(new_dict[11]), '..')


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
