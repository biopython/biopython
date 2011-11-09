# Copyright 2010-2011 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Test code for working with BGZF files (used in BAM files).

See also the doctests in bgzf.py which are called via run_tests.py
"""

import unittest
import gzip
import os

from Bio._py3k import _as_bytes
_empty_bytes_string = _as_bytes("")

from Bio import bgzf

class BgzfTests(unittest.TestCase):
    def setUp(self):
        self.temp_file = "temp.bgzf"
        if os.path.isfile(self.temp_file):
            os.remove(self.temp_file)

    def tearDown(self):
        if os.path.isfile(self.temp_file):
            os.remove(self.temp_file)

    def rewrite(self, compressed_input_file, output_file):
        h = gzip.open(compressed_input_file, "rb")
        data = h.read()
        h.close()

        h = bgzf.BgzfWriter(output_file, "wb")
        h.write(data)
        h.flush()
        h.flush() #Second flush gives empty BGZF block as BAM EOF marker
        h.close()

        h = gzip.open(output_file)
        new_data = h.read()
        h.close()

        #Check the decompressed files agree
        self.assert_(new_data, "Empty BGZF file?")
        self.assertEqual(len(data), len(new_data))
        self.assertEqual(data, new_data)

    def check_blocks(self, old_file, new_file):
        h = open(old_file, "rb")
        old = list(bgzf.BgzfBlocks(h))
        h.close()
        h = open(new_file, "rb")
        new = list(bgzf.BgzfBlocks(h))
        h.close()
        self.assertEqual(len(old), len(new))
        self.assertEqual(old, new)

    def check_random(self, filename):
        """Check BGZF random access by reading blocks in forward & reverse order"""
        h = gzip.open(filename, "rb")
        old = h.read()
        h.close()

        h = open(filename, "rb")
        blocks = list(bgzf.BgzfBlocks(h))
        h.close()

        #Forward
        new = _empty_bytes_string
        h = bgzf.BgzfReader(filename, "rb")
        for start, raw_len, data_len in blocks:
            #print start, raw_len, data_len
            h.seek(bgzf.make_virtual_offset(start,0))
            data = h.read(data_len)
            #self.assertEqual(start + raw_len, h._handle.tell())
            new += data
        h.close()
        self.assertEqual(len(old), len(new))
        self.assertEqual(old, new)

        #Reverse
        new = _empty_bytes_string
        h = bgzf.BgzfReader(filename, "rb")
        for start, raw_len, data_len in blocks[::-1]:
            #print start, raw_len, data_len
            h.seek(bgzf.make_virtual_offset(start,0))
            data = h.read(data_len)
            #self.assertEqual(start + raw_len, h._handle.tell())
            new = data + new
        h.close()
        self.assertEqual(len(old), len(new))
        self.assertEqual(old, new)



    def test_random_bam_ex1(self):
        """Check random access to SamBam/ex1.bam"""
        self.check_random("SamBam/ex1.bam")

    def test_random_bam_ex1_refresh(self):
        """Check random access to SamBam/ex1_refresh.bam"""
        self.check_random("SamBam/ex1_refresh.bam")

    def test_random_bam_ex1_header(self):
        """Check random access to SamBam/ex1_header.bam"""
        self.check_random("SamBam/ex1_header.bam")

    def test_random_example_fastq(self):
        """Check random access to Quality/example.fastq.bgz"""
        self.check_random("Quality/example.fastq.bgz")

    def test_bam_ex1(self):
        """Reproduce BGZF compression for BAM file"""
        temp_file = self.temp_file

        #Note this example is from an old version of samtools
        #and all the blocks are full (except the last one)
        self.rewrite("SamBam/ex1.bam", temp_file)

        #Now check the blocks agree (using the fact that
        #this example BAM file has simple block usage)
        self.check_blocks("SamBam/ex1.bam", temp_file)

    def test_example_fastq(self):
        """Reproduce BGZF compression for a FASTQ file"""
        temp_file = self.temp_file
        self.rewrite("Quality/example.fastq.gz", temp_file)
        self.check_blocks("Quality/example.fastq.bgz", temp_file)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
