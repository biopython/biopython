# Copyright 2006 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
import unittest
from Bio.PopGen import GenePop
from Bio.PopGen import FDist
from Bio.PopGen.GenePop import FileParser
from Bio.PopGen.FDist.Utils import convert_genepop_to_fdist

#Tests fdist related code. Note: this case doesn't require fdist
#test_PopGen_FDist tests code that requires fdist


class RecordTest(unittest.TestCase):
    def test_record_basic(self):
        """Basic test on Record
        """

        r = FDist.Record()
        assert type(r.data_org)  == int
        assert type(r.num_pops)  == int
        assert type(r.num_loci)  == int
        assert type(r.loci_data) == list

class ParserTest(unittest.TestCase):
    def setUp(self):
        files = ["fdist1"]
        self.handles = []
        for filename in files:
            self.handles.append(open(os.path.join("PopGen", filename)))

        self.pops_loci = [
            (3, 4)
        ]
        self.num_markers = [
            [2, 3, 4, 2]
        ]
        #format is locus, pop, position, value
        self.test_pos = [
            [
                (0, 0, 0, 5),
                (3, 2, 0, 5)
            ]
        ]

    def tearDown(self):
        for handle in self.handles:
            handle.close()

    def test_record_parser(self):
        """Basic operation of the Record Parser.
        """
        for index in range(len(self.handles)):
            handle = self.handles[index]
            rec = FDist.read(handle)
            assert isinstance(rec, FDist.Record)
            assert rec.data_org == 0 #We don't support any other
            assert rec.num_pops, rec.num_loci == self.pops_loci[index]
            for i in range(len(self.num_markers[index])):
                assert rec.loci_data[i][0] == \
                       self.num_markers[index][i]
            for i in range(len(self.test_pos[index])):
                my_test_pos = self.test_pos[index]
                for test in my_test_pos:
                    locus, pop, pos, value = test
                    assert(rec.loci_data[locus][1][pop][pos] == value)

class ConversionTest(unittest.TestCase):
    def setUp(self):
        files = ["c2line.gen", "haplo2.gen"]
        self.handles = []
        self.names = []
        for filename in files:
            self.names.append(os.path.join("PopGen", filename))
            self.handles.append(open(os.path.join("PopGen", filename)))

    def test_convert(self):
        """Basic conversion test.
        """
        for i in range(len(self.names)):
            handle = self.handles[i]
            gp_rec = GenePop.read(handle)
            fd_rec = convert_genepop_to_fdist(gp_rec)
            assert(fd_rec.num_loci == 3)
            assert(fd_rec.num_pops == 3)

    def test_convert_big(self):
        """Big interface conversion test.
        """
        for i in range(len(self.names)):
            gp_rec = FileParser.read(self.names[i])
            fd_rec = convert_genepop_to_fdist(gp_rec)
            assert(fd_rec.num_loci == 3)
            assert(fd_rec.num_pops == 3)
            gp_rec._handle.close() #TODO - Needs a proper fix

    def tearDown(self):
        for handle in self.handles:
            handle.close()

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
