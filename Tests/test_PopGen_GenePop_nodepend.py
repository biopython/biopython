# Copyright 2007 by Tiago Antao <tiagoantao@gmail.com>.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


import os
import unittest
from Bio.PopGen import GenePop
from Bio.PopGen.GenePop import FileParser


class RecordTest(unittest.TestCase):
    def test_record_basic(self):
        """Basic test on Record
        """

        r = GenePop.Record()
        assert type(r.marker_len)   == int
        assert type(r.comment_line) == str
        assert type(r.loci_list)    == list
        assert type(r.populations)  == list


class ParserTest(unittest.TestCase):
    def setUp(self):
        files = ["c2line.gen", "c3line.gen", "c2space.gen", "c3space.gen",
                 "haplo3.gen", "haplo2.gen"]
        self.handles = []
        for filename in files:
            self.handles.append(open(os.path.join("PopGen", filename)))

        self.pops_indivs = [
            (3, [4, 3, 5]),
            (3, [4, 3, 5]),
            (3, [4, 3, 5]),
            (3, [4, 3, 5]),
            (3, [4, 3, 5]),
            (3, [4, 3, 5])
        ]
        self.num_loci = [3, 3, 3, 3, 3, 3]
        self.marker_len = [2, 3, 2, 3, 3, 2]
        self.pop_names = ["4", "b3", "5"]

    def tearDown(self):
        for handle in self.handles:
            handle.close()

    def test_record_parser(self):
        """Basic operation of the Record Parser.
        """
        for index in range(len(self.handles)):
            handle = self.handles[index]
            rec = GenePop.read(handle)
            assert isinstance(rec, GenePop.Record)
            assert len(rec.loci_list) == self.num_loci[index]
            assert rec.marker_len == self.marker_len[index]
            assert len(rec.populations) == self.pops_indivs[index][0]
            assert rec.pop_list == self.pop_names
            for i in range(self.pops_indivs[index][0]):
                assert len(rec.populations[i]) == \
                           self.pops_indivs[index][1][i]

    def test_wrong_file_parser(self):
        """Testing the ability to deal with wrongly formatted files
        """
        f = open(os.path.join("PopGen", "fdist1"))
        try:
            rec = GenePop.read(f)
            raise Error("Should have raised exception")
        except ValueError:
            pass
        f.close()


class FileParserTest(unittest.TestCase):
    def setUp(self):
        self.files = map(lambda x: os.path.join("PopGen", x),
             ["c2line.gen", "c3line.gen", "c2space.gen",
              "c3space.gen", "haplo3.gen", "haplo2.gen"])
        self.pops_indivs = [
            (3, [4, 3, 5]),
            (3, [4, 3, 5]),
            (3, [4, 3, 5]),
            (3, [4, 3, 5]),
            (3, [4, 3, 5]),
            (3, [4, 3, 5])
        ]
        self.num_loci = [3, 3, 3, 3, 3, 3]

    def test_file_record_parser(self):
        """Basic operation of the File Record Parser.
        """
        for index in range(len(self.files)):
            fname = self.files[index]
            rec = FileParser.read(fname)
            assert isinstance(rec, FileParser.FileRecord)
            assert len(rec.loci_list) == self.num_loci[index]
            for skip in range(self.pops_indivs[index][0]):
                if rec.skip_population() == False:
                    raise Error("Not enough populations")
            if rec.skip_population() == True:
                    raise Error("Too much populations")
            for i in range(self.pops_indivs[index][0]):
                continue
                assert len(rec.populations[i]) == \
                           self.pops_indivs[index][1][i]
            rec._handle.close() #TODO - Needs a proper fix

    def test_wrong_file_parser(self):
        """Testing the ability to deal with wrongly formatted files
        """
        f = open(os.path.join("PopGen", "fdist1"))
        try:
            rec = GenePop.read(f)
            raise Error("Should have raised exception")
        except ValueError:
            pass
        f.close()

class UtilsTest(unittest.TestCase):
    def setUp(self):
        #All files have to have at least 3 loci and 2 pops
        files = ["c2line.gen"]
        self.handles = []
        for filename in files:
            self.handles.append(open(os.path.join("PopGen", filename)))


    def tearDown(self):
        for handle in self.handles:
            handle.close()

    def test_utils(self):
        """Basic operation of GenePop Utils.
        """
        for index in range(len(self.handles)):
            handle = self.handles[index]
            rec = GenePop.read(handle)
        initial_pops = len(rec.populations)
        initial_loci = len(rec.loci_list)
        first_loci = rec.loci_list[0]
        rec.remove_population(0)
        assert len(rec.populations) == initial_pops - 1
        rec.remove_locus_by_name(first_loci)
        assert len(rec.loci_list) == initial_loci - 1
        assert rec.loci_list[0] != first_loci
        rec.remove_locus_by_position(0)
        assert len(rec.loci_list) == initial_loci - 2

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
