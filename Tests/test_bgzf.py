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
from random import shuffle

from Bio._py3k import _as_bytes, _as_string
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
        self.assertFalse(h.seekable())
        self.assertFalse(h.isatty())
        self.assertEqual(h.fileno(), h._handle.fileno())
        h.close()  # Gives empty BGZF block as BAM EOF marker

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

    def check_text(self, old_file, new_file):
        h = open(old_file)  # text mode!
        old_line = h.readline()
        old = old_line + h.read()
        h.close()

        h = bgzf.BgzfReader(new_file, "r")  # Text mode!
        new_line = h.readline()
        new = new_line + h.read(len(old))
        h.close()

        self.assertEqual(old_line, new_line)
        self.assertEqual(len(old), len(new))
        self.assertEqual(old, new)

    def check_by_line(self, old_file, new_file, old_gzip=False):
        for mode in ["r", "rb"]:
            if old_gzip:
                h = gzip.open(old_file, mode)
            else:
                h = open(old_file, mode)
            old = h.read()
            #Seems gzip can return bytes even if mode="r",
            #perhaps a bug in Python 3.2?
            if "b" in mode:
                old = _as_bytes(old)
            else:
                old = _as_string(old)
            h.close()

            for cache in [1,10]:
                h = bgzf.BgzfReader(new_file, mode, max_cache=cache)
                if "b" in mode:
                    new = _empty_bytes_string.join(line for line in h)
                else:
                    new = "".join(line for line in h)
                h.close()

                self.assertEqual(len(old), len(new))
                self.assertEqual(old[:10], new[:10],
                                 "%r vs %r, mode %r" % (old[:10], new[:10], mode))
                self.assertEqual(old, new)

    def check_by_char(self, old_file, new_file, old_gzip=False):
        for mode in ["r", "rb"]:
            if old_gzip:
                h = gzip.open(old_file,mode)
            else:
                h = open(old_file, mode)
            old = h.read()
            #Seems gzip can return bytes even if mode="r",
            #perhaps a bug in Python 3.2?
            if "b" in mode:
                old = _as_bytes(old)
            else:
                old = _as_string(old)
            h.close()

            for cache in [1,10]:
                h = bgzf.BgzfReader(new_file, mode, max_cache=cache)
                temp = []
                while True:
                    char = h.read(1)
                    if not char:
                        break
                    temp.append(char)
                if "b" in mode:
                    new = _empty_bytes_string.join(temp)
                else:
                    new = "".join(temp)
                del temp
                h.close()

                self.assertEqual(len(old), len(new))
                #If bytes vs unicode mismatch, give a short error message:
                self.assertEqual(old[:10], new[:10],
                                 "%r vs %r, mode %r" % (old[:10], new[:10], mode))
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
        self.assertTrue(h.seekable())
        self.assertFalse(h.isatty())
        self.assertEqual(h.fileno(), h._handle.fileno())
        for start, raw_len, data_start, data_len in blocks:
            #print start, raw_len, data_start, data_len
            h.seek(bgzf.make_virtual_offset(start,0))
            data = h.read(data_len)
            self.assertEqual(len(data), data_len)
            #self.assertEqual(start + raw_len, h._handle.tell())
            self.assertEqual(len(new), data_start)
            new += data
        h.close()
        self.assertEqual(len(old), len(new))
        self.assertEqual(old, new)

        #Reverse
        new = _empty_bytes_string
        h = bgzf.BgzfReader(filename, "rb")
        for start, raw_len, data_start, data_len in blocks[::-1]:
            #print start, raw_len, data_start, data_len
            h.seek(bgzf.make_virtual_offset(start,0))
            data = h.read(data_len)
            self.assertEqual(len(data), data_len)
            #self.assertEqual(start + raw_len, h._handle.tell())
            new = data + new
        h.close()
        self.assertEqual(len(old), len(new))
        self.assertEqual(old, new)

        #Jump back - non-sequential seeking
        if len(blocks) >= 3:
            h = bgzf.BgzfReader(filename, "rb", max_cache = 1)
            #Seek to a late block in the file,
            #half way into the third last block
            start, raw_len, data_start, data_len = blocks[-3]
            voffset = bgzf.make_virtual_offset(start, data_len // 2)
            h.seek(voffset)
            self.assertEqual(voffset, h.tell())
            data = h.read(1000)
            self.assertTrue(data in old)
            self.assertEqual(old.find(data), data_start + data_len // 2)
            #Now seek to an early block in the file,
            #half way into the second block
            start, raw_len, data_start, data_len = blocks[1]
            h.seek(bgzf.make_virtual_offset(start, data_len // 2))
            voffset = bgzf.make_virtual_offset(start, data_len // 2)
            h.seek(voffset)
            self.assertEqual(voffset, h.tell())
            #Now read all rest of this block and start of next block
            data = h.read(data_len + 1000)
            self.assertTrue(data in old)
            self.assertEqual(old.find(data), data_start + data_len // 2)
            h.close()

        #Check seek/tell at block boundaries
        v_offsets = []
        for start, raw_len, data_start, data_len in blocks:
            for within_offset in [0, 1, data_len // 2, data_len - 1]:
                if within_offset < 0 or data_len <= within_offset:
                    continue
                voffset = bgzf.make_virtual_offset(start, within_offset)
                real_offset = data_start + within_offset
                v_offsets.append((voffset, real_offset))
        shuffle(v_offsets)
        h = bgzf.BgzfReader(filename, "rb", max_cache = 1)
        for voffset, real_offset in v_offsets:
            h.seek(0)
            self.assertTrue(voffset >= 0 and real_offset >= 0)
            self.assertEqual(h.read(real_offset), old[:real_offset])
            self.assertEqual(h.tell(), voffset)
        for voffset, real_offset in v_offsets:
            h.seek(voffset)
            self.assertEqual(h.tell(), voffset)
        h.close()

    def test_random_bam_ex1(self):
        """Check random access to SamBam/ex1.bam"""
        self.check_random("SamBam/ex1.bam")

    def test_random_bam_ex1_refresh(self):
        """Check random access to SamBam/ex1_refresh.bam"""
        self.check_random("SamBam/ex1_refresh.bam")

    def test_random_bam_ex1_header(self):
        """Check random access to SamBam/ex1_header.bam"""
        self.check_random("SamBam/ex1_header.bam")

    def test_random_wnts_xml(self):
        """Check random access to Blast/wnts.xml.bgz"""
        self.check_random("Blast/wnts.xml.bgz")

    def test_random_example_fastq(self):
        """Check random access to Quality/example.fastq.bgz"""
        self.check_random("Quality/example.fastq.bgz")

    def test_random_example_cor6(self):
        """Check random access to GenBank/cor6_6.gb.bgz"""
        self.check_random("GenBank/cor6_6.gb.bgz")

    def test_text_wnts_xml(self):
        """Check text mode access to Blast/wnts.xml.bgz"""
        self.check_text("Blast/wnts.xml", "Blast/wnts.xml.bgz")

    def test_text_example_fastq(self):
        """Check text mode access to Quality/example.fastq.bgz"""
        self.check_text("Quality/example.fastq", "Quality/example.fastq.bgz")

    def test_iter_wnts_xml(self):
        """Check iteration over Blast/wnts.xml.bgz"""
        self.check_by_line("Blast/wnts.xml", "Blast/wnts.xml.bgz")
        self.check_by_char("Blast/wnts.xml", "Blast/wnts.xml.bgz")

    def test_iter_example_fastq(self):
        """Check iteration over Quality/example.fastq.bgz"""
        self.check_by_line("Quality/example.fastq", "Quality/example.fastq.bgz")
        self.check_by_char("Quality/example.fastq", "Quality/example.fastq.bgz")

    def test_iter_example_cor6(self):
        """Check iteration over GenBank/cor6_6.gb.bgz"""
        self.check_by_line("GenBank/cor6_6.gb", "GenBank/cor6_6.gb.bgz")
        self.check_by_char("GenBank/cor6_6.gb", "GenBank/cor6_6.gb.bgz")

    def test_iter_example_gb(self):
        """Check iteration over GenBank/NC_000932.gb.bgz"""
        self.check_by_line("GenBank/NC_000932.gb", "GenBank/NC_000932.gb.bgz")
        self.check_by_char("GenBank/NC_000932.gb", "GenBank/NC_000932.gb.bgz")

    def test_bam_ex1(self):
        """Reproduce BGZF compression for BAM file"""
        temp_file = self.temp_file

        #Note this example is from an old version of samtools
        #and all the blocks are full (except the last one)
        self.rewrite("SamBam/ex1.bam", temp_file)

        #Now check the blocks agree (using the fact that
        #this example BAM file has simple block usage)
        self.check_blocks("SamBam/ex1.bam", temp_file)

    def test_iter_bam_ex1(self):
        """Check iteration over SamBam/ex1.bam"""
        self.check_by_char("SamBam/ex1.bam", "SamBam/ex1.bam", True)

    def test_example_fastq(self):
        """Reproduce BGZF compression for a FASTQ file"""
        temp_file = self.temp_file
        self.rewrite("Quality/example.fastq.gz", temp_file)
        self.check_blocks("Quality/example.fastq.bgz", temp_file)

    def test_example_gb(self):
        """Reproduce BGZF compression for NC_000932 GenBank file"""
        temp_file = self.temp_file
        self.rewrite("GenBank/NC_000932.gb.bgz", temp_file)
        self.check_blocks("GenBank/NC_000932.gb.bgz", temp_file)

    def test_example_cor6(self):
        """Reproduce BGZF compression for cor6_6.gb GenBank file"""
        temp_file = self.temp_file
        self.rewrite("GenBank/cor6_6.gb.bgz", temp_file)
        self.check_blocks("GenBank/cor6_6.gb.bgz", temp_file)

    def test_example_wnts_xml(self):
        """Reproduce BGZF compression for wnts.xml BLAST file"""
        temp_file = self.temp_file
        self.rewrite("Blast/wnts.xml.bgz", temp_file)
        self.check_blocks("Blast/wnts.xml.bgz", temp_file)

    def test_write_tell(self):
        """Check offset works during BGZF writing"""
        temp_file = self.temp_file

        h = bgzf.open(temp_file, "w") #Text mode!
        h.write("X" * 100000)
        offset = h.tell()
        self.assertNotEqual(offset, 100000) #Should be a virtual offset!
        h.write("Magic" + "Y" * 100000)
        h.close()

        h = bgzf.open(temp_file, "r") #Text mode!
        h.seek(offset)
        self.assertEqual(h.read(5), "Magic")
        h.close()


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
