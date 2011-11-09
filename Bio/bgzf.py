#!/usr/bin/env python
# Copyright 2010-2011 by Peter Cock.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
r"""Fairly low level API for working with BGZF files (e.g. BAM files).

The SAM/BAM file format (Sequence Alignment/Map) comes in a plain text
format (SAM), and a compressed binary format (BAM). The latter uses a
modified form of gzip compression called BGZF, which in principle can
be applied to any file format. This is described together with the
SAM/BAM file format at http://samtools.sourceforge.net/SAM1.pdf


Aim of this module
------------------

The Python gzip library can be used to read BGZF files, since for
decompression they are just (specialised) gzip files. What this
module aims to facilitate is random access to BGZF files (using the
'virtual offset' idea), and writing BGZF files (which means using
suitably sized gzip blocks and writing the extra 'BC' field in the
gzip headers). The existing gzip library will be used internally.

Initially this will be used to provide random access and writing of
BAM files. However, the BGZF format could also be used on other
sequential data (in the sense of one record after another), such
as most of the sequence data formats supported in Bio.SeqIO (like
FASTA, FASTQ, GenBank, etc).


Technical Introduction to BGZF
------------------------------

The gzip file format allows multiple compressed blocks, each of which
could be a stand alone gzip file. As an interesting bonus, this means
you can use Unix "cat" to combined to gzip files into one by
concatenating them. Also, each block can have one of several compression
levels (including uncompressed, which actually takes up a little bit
more space due to the gzip header).

What the BAM designers realised was that random access to data stored
in traditional gzip files was slow, breaking the file into gzip blocks
would allow fast random access to each block. To access a particular
piece of the decompressed data, you just need to know which block it
starts in (the offset of the gzip block start), and how far into the
(decompressed) contents of the block you need to read.

One problem with this is finding the gzip block sizes efficiently.
You can do it with a standard gzip file, but it requires every block
to be decompressed -- and that would be rather slow.

All that differs in BGZF is that compressed size of each gzip block
is limited to 2^16 bytes, and an extra 'BC' field in the gzip header
records this size. Traditional decompression tools can ignore this,
and unzip the file just like any other gzip file.

The point of this is you can look at the first BGZF block, find out
how big it is from this 'BC' header, and thus seek immediately to
the second block, and so on.

The BAM indexing scheme records read positions using a 64 bit
'virtual offset', comprising coffset<<16|uoffset, where uoffset is
the file offset of the BGZF block containing the start of the read
(unsigned integer using up to 64-16 = 48 bits), and coffset is the
offset within the (decompressed) block (unsigned 16 bit integer).

This limits you to BAM files where the last block starts by 2^48
bytes, or 256 petabytes, and the decompressed size of each block
is at most 2^16 bytes, or 64kb. Note that this matches the BGZF
'BC' field size which limits the compressed size of each block to
2^16 bytes, allowing for BAM files to use BGZF with no gzip
compression (useful for intermediate files in memory to reduced
CPU load).


Example
-------

This is an ordinary GenBank file compressed using BGZF, so it can
be decompressed using gzip,

>>> import gzip
>>> handle = gzip.open("GenBank/NC_000932.gb.bgz", "rb")
>>> handle.tell()
0
>>> handle.readline()
'LOCUS       NC_000932             154478 bp    DNA     circular PLN 15-APR-2009\n'
>>> handle.tell()
80
>>> handle.readline()
'DEFINITION  Arabidopsis thaliana chloroplast, complete genome.\n'
>>> handle.tell()
143
>>> data = handle.read(70000)
>>> handle.tell()
70143
>>> handle.readline()
'f="GeneID:844718"\n'
>>> handle.readline()
'     CDS             complement(join(84337..84771,85454..85843))\n'
>>> offset = handle.seek(65536*3 + 126)
>>> handle.readline()
'    68521 tatgtcattc gaaattgtat aaagacaact cctatttaat agagctattt gtgcaagtat\n'
>>> handle.close()

We can also access the file using the BGZF reader - but pay
attention to the file offsets which will be explained below:

>>> handle = BgzfReader("GenBank/NC_000932.gb.bgz", "rb")
>>> handle.tell()
0
>>> handle.readline()
'LOCUS       NC_000932             154478 bp    DNA     circular PLN 15-APR-2009\n'
>>> handle.tell()
80
>>> handle.readline()
'DEFINITION  Arabidopsis thaliana chloroplast, complete genome.\n'
>>> handle.tell()
143
>>> data = handle.read(70000)
>>> handle.tell()
987828735
>>> handle.readline()
'f="GeneID:844718"\n'
>>> handle.readline()
'     CDS             complement(join(84337..84771,85454..85843))\n'
>>> offset = handle.seek(make_virtual_offset(55074, 126))
>>> handle.readline()
'    68521 tatgtcattc gaaattgtat aaagacaact cctatttaat agagctattt gtgcaagtat\n'
>>> handle.close()

Notice the handle's offset looks different as a BGZF file. This
brings us to the key point about BGZF, which is the block structure:

>>> handle = open("GenBank/NC_000932.gb.bgz", "rb")
>>> for values in BgzfBlocks(handle):
...     print "Start %i, length %i, data %i" % values
Start 0, length 15073, data 65536
Start 15073, length 17857, data 65536
Start 32930, length 22144, data 65536
Start 55074, length 22230, data 65536
Start 77304, length 14939, data 43478
Start 92243, length 28, data 0
>>> handle.close()                             

By reading ahead 70,000 bytes we moved into the second BGZF block,
and at that point the BGZF virtual offsets start to look different
a simple offset into the decompressed data as exposed by the gzip
library.

Using the seek for the decompressed co-ordinates, 65536*3 + 126
is equivalent to jumping the first thre blocks (each size 65536
after decompression) and starting at byte 126 of the third block
(after decompression). For BGZF, we need to know the block's
offset of 55074 and the offset within the block of 126 to get
the BGZF virtual offset.

The catch with BGZF virtual offsets is while they can be compared
(which offset comes first in the file), you cannot safely subtract
them to get the size of the data between them, nor add/subtract
a relative offset.

Of course you can parse this file with Bio.SeqIO using BgzfReader,
although there isn't any benefit over using gzip.open(...), unless
you want to index BGZF compressed sequence files:

>>> from Bio import SeqIO
>>> handle = BgzfReader("GenBank/NC_000932.gb.bgz", "rb")
>>> record = SeqIO.read(handle, "genbank")
>>> handle.close()
>>> print record.id
NC_000932.1

"""
#TODO - Move somewhere else in Bio.* namespace?

import gzip
import zlib
import struct
import __builtin__ #to access the usual open function

def bgzf_open(filename, mode):
    if "r" in mode.lower():
        return BgzfReader(filename, mode)
    elif "w" in mode.lower() or "a" in mode.lower():
        return BgzfWriter(filename, mode)
    else:
        raise ValueError("Bad mode %r" % mode)

def make_virtual_offset(block_start_offset, within_block_offset):
    """Compute a BGZF virtual offset from block start and within block offsets.

    The BAM indexing scheme records read positions using a 64 bit
    'virtual offset', comprising in C terms:

    within_block_offset<<16 | block_start_offset

    Here block_start_offset is the file offset of the BGZF block
    start (unsigned integer using up to 64-16 = 48 bits), and
    within_block_offset within the (decompressed) block (unsigned
    16 bit integer).

    >>> make_virtual_offset(0,0)
    0
    >>> make_virtual_offset(0,1)
    1
    >>> make_virtual_offset(0, 2**16 - 1)
    65535
    >>> make_virtual_offset(0, 2**16)
    Traceback (most recent call last):
    ...
    ValueError: Require 0 <= within_block_offset < 2**16, got 65536

    >>> make_virtual_offset(1,0)
    65536
    >>> make_virtual_offset(1,1)
    65537
    >>> make_virtual_offset(1, 2**16 - 1)
    131071

    >>> make_virtual_offset(100000,0)
    6553600000
    >>> make_virtual_offset(100000,1)
    6553600001
    >>> make_virtual_offset(100000,10)
    6553600010

    >>> make_virtual_offset(2**48,0)
    Traceback (most recent call last):
    ...
    ValueError: Require 0 <= block_start_offset < 2**48, got 281474976710656

    """
    if within_block_offset < 0 or within_block_offset >= 2**16:
        raise ValueError("Require 0 <= within_block_offset < 2**16, got %r" % within_block_offset)
    if block_start_offset < 0 or block_start_offset >= 2**48:
        raise ValueError("Require 0 <= block_start_offset < 2**48, got %r" % block_start_offset)
    return (block_start_offset<<16) | within_block_offset

def split_virtual_offset(virtual_offset):
    """Divides a 64-bit BGZF virtual offset into block start & within block offsets.

    >>> split_virtual_offset(6553600000)
    (100000, 0)
    >>> split_virtual_offset(6553600010)
    (100000, 10)

    """
    start = virtual_offset>>16
    return start, virtual_offset ^ (start<<16)

def BgzfBlocks(handle):
    """Low level debugging function to inspect BGZF blocks.

    Returns the block start offset (see virtual offsets), the block
    length (add these for the start of the next block), and the
    decompressed length of the blocks contents (limited to 65536 in
    BGZF).

    >>> handle = open("SamBam/ex1.bam", "rb")
    >>> for values in BgzfBlocks(handle):
    ...     print "Start %i, length %i, data %i" % values
    Start 0, length 18239, data 65536
    Start 18239, length 18223, data 65536
    Start 36462, length 18017, data 65536
    Start 54479, length 17342, data 65536
    Start 71821, length 17715, data 65536
    Start 89536, length 17728, data 65536
    Start 107264, length 17292, data 63398
    Start 124556, length 28, data 0
    >>> handle.close()

    Indirectly we can tell this file came from an old version of
    samtools because all the blocks (except the final one and the
    dummy empty EOF marker block) are 65536 bytes.  Later versions
    avoid splitting a read between two blocks, and give the header
    its own block (useful to speed up replacing the header). You
    can see this in ex1_refresh.bam created using samtools 0.1.18:

    samtools view -b ex1.bam > ex1_refresh.bam

    >>> handle = open("SamBam/ex1_refresh.bam", "rb")
    >>> for values in BgzfBlocks(handle):
    ...     print "Start %i, length %i, data %i" % values
    Start 0, length 53, data 38
    Start 53, length 18195, data 65434
    Start 18248, length 18190, data 65409
    Start 36438, length 18004, data 65483
    Start 54442, length 17353, data 65519
    Start 71795, length 17708, data 65411
    Start 89503, length 17709, data 65466
    Start 107212, length 17390, data 63854
    Start 124602, length 28, data 0
    >>> handle.close()

    The above example has no embedded SAM header (thus the first block
    is very small), while the next example does. Notice that the rest
    of the blocks show the same sizes (the contain the same read data):

    >>> handle = open("SamBam/ex1_header.bam", "rb")
    >>> for values in BgzfBlocks(handle):
    ...     print "Start %i, length %i, data %i" % values
    Start 0, length 104, data 103
    Start 104, length 18195, data 65434
    Start 18299, length 18190, data 65409
    Start 36489, length 18004, data 65483
    Start 54493, length 17353, data 65519
    Start 71846, length 17708, data 65411
    Start 89554, length 17709, data 65466
    Start 107263, length 17390, data 63854
    Start 124653, length 28, data 0
    >>> handle.close()

    """
    while True:
        start_offset = handle.tell()
        #This may raise StopIteration which is perfect here
        block_length, data = _load_bgzf_block(handle)
        yield start_offset, block_length, len(data)


def _load_bgzf_block(handle):
    #Change indentation later...
    magic = handle.read(4)
    if not magic:
        #End of file
        raise StopIteration
    if magic != "\x1f\x8b\x08\x04":
        raise ValueError(r"A BGZF (e.g. a BAM file) block should start with "
                         r"'\x1f\x8b\x08\x04' (decimal 31 139 8 4), not %s."
                         r"handle.tell() now says %r"
                         % (repr(magic), handle.tell()))
    gzip_mod_time = handle.read(4) #uint32_t
    gzip_extra_flags = handle.read(1) #uint8_t
    gzip_os = handle.read(1) #uint8_t
    extra_len = struct.unpack("<H", handle.read(2))[0] #uint16_t
        
    block_size = None
    x_len = 0
    while x_len < extra_len:
        subfield_id = handle.read(2)
        subfield_len = struct.unpack("<H", handle.read(2))[0] #uint16_t
        subfield_data = handle.read(subfield_len)
        x_len += subfield_len + 4
        if subfield_id == "BC":
            assert subfield_len == 2, "Wrong BC payload length"
            assert block_size is None, "Two BC subfields?"
            block_size = struct.unpack("<H", subfield_data)[0]+1 #uint16_t
    assert x_len == extra_len, (x_len, extra_len)
    assert block_size is not None, "Missing BC, this isn't a BGZF file!"
    #Now comes the compressed data, CRC, and length of uncompressed data.
    deflate_size = block_size - 1 - extra_len - 19
    d = zlib.decompressobj(-15) #Negative window size means no headers
    data = d.decompress(handle.read(deflate_size)) + d.flush()
    expected_crc = handle.read(4)
    expected_size = struct.unpack("<I", handle.read(4))[0]
    assert expected_size == len(data), \
           "Decompressed to %i, not %i" % (len(data), expected_size)
    #Should cope with a mix of Python platforms...
    crc = zlib.crc32(data)
    if crc < 0:
        crc = struct.pack("<i", crc)
    else:
        crc = struct.pack("<I", crc)
    assert expected_crc == crc, \
           "CRC is %s, not %s" % (crc, expected_crc)
    return block_size, data


class BgzfReader(object):
    r"""BGZF reader, acts like a read only handle but seek/tell differ.

    Let's use the BgzfBlocks function to have a peak at the BGZF blocks
    in an example BAM file,

    >>> handle = open("SamBam/ex1.bam", "rb")
    >>> for values in BgzfBlocks(handle):
    ...     print "Start %i, length %i, data %i" % values
    Start 0, length 18239, data 65536
    Start 18239, length 18223, data 65536
    Start 36462, length 18017, data 65536
    Start 54479, length 17342, data 65536
    Start 71821, length 17715, data 65536
    Start 89536, length 17728, data 65536
    Start 107264, length 17292, data 63398
    Start 124556, length 28, data 0
    >>> handle.close()
 
    Now let's see how to use this block information to jump to
    specific parts of the decompressed BAM file:

    >>> handle = BgzfReader("SamBam/ex1.bam", "rb")
    >>> handle.tell()
    0
    >>> handle.read(4)
    'BAM\x01'
    >>> handle.tell()
    4

    So far nothing so strange, we got the magic marker used at the
    start of a decompressed BAM file, and the handle position makes
    sense. Now however, let's jump to the end of this block and 4
    bytes into the next block by reading 65536 bytes,

    >>> data = handle.read(65536)
    >>> len(data)
    65536
    >>> handle.tell()
    1195311108

    Expecting 4 + 65536 = 65540 were you? Well this is a BGZF 64-bit
    virtual offset, which means:

    >>> split_virtual_offset(1195311108)
    (18239, 4)

    You should spot 18239 as the start of the second BGZF block, while
    the 4 is the offset into this block. See also make_virtual_offset,

    >>> make_virtual_offset(18239, 4)
    1195311108

    Let's jump back to almost the start of the file,

    >>> make_virtual_offset(0, 2)
    2
    >>> handle.seek(2)
    2
    >>> handle.read(2)
    'M\x01'
    >>> handle.close()

    Note that you can use the max_cache argument to limit the number of
    BGZF blocks cached in memory. The default is 100, and since each
    block can be up to 64kb, the default cache could take up to 6MB of
    RAM. The cache is not important for reading through the file in one
    pass, but is important for improving performance of random access.
    """

    def __init__(self, filename=None, mode=None, fileobj=None, max_cache=100):
        if fileobj:
            assert filename is None and mode is None
            handle = fileobj
        else:
            if "w" in mode.lower() \
            or "a" in mode.lower():
                raise ValueError("Must use read mode (default), not write or append mode")
            handle = __builtin__.open(filename, mode)
        self._handle = handle
        self.max_cache = max_cache
        self._buffers = {}
        self._block_start_offset = None
        self._load_block()

    def _load_block(self, start_offset=None):
        if start_offset is None:
            start_offset = self._handle.tell()
        if start_offset == self._block_start_offset:
            self._within_block_offset = 0
            return
        elif start_offset in self._buffers:
            #Already in cache
            self._block_size, self._buffer = self._buffers[start_offset]
            self._within_block_offset = 0
            return
        #Must hit the disk... first check cache limits,
        while len(self._buffers) >= self.max_cache:
            #TODO - Implemente LRU cache removal?
            self._buffers.pop(0)
        #Now load the block
        handle = self._handle
        if start_offset is not None:
            handle.seek(start_offset)
        self._block_start_offset = handle.tell()
        try:
            self._block_size, self._buffer = _load_bgzf_block(handle)
        except StopIteration:
            #EOF
            self._block_size = 0
            self._buffer = ""
        self._within_block_offset = 0
        #Finally save the block in our cache,
        self._buffers[self._block_start_offset] = (self._block_size, self._buffer)

    def tell(self):
        """Returns a 64-bit unsigned BGZF virtual offset."""
        if self._within_block_offset == 65536:
            assert 65536 == len(self._buffer)
            #At the end of one block, start of next
            #TODO - Handle this in the read method?
            self._load_block()
        return make_virtual_offset(self._block_start_offset, self._within_block_offset)

    def seek(self, virtual_offset):
        """Seek to a 64-bit unsigned BGZF virtual offset."""
        start_offset, within_block = split_virtual_offset(virtual_offset)
        if start_offset != self._block_start_offset:
            #Don't need to load the block if already there
            #(this avoids a function call since _load_block would do nothing)
            self._load_block(start_offset)
        if within_block >= len(self._buffer) \
        and not (within_block == 0 and len(self._buffer)==0):
            raise ValueError("Within offset %i but block size only %i" \
                             % (within_block, len(self._buffer)))
        self._within_block_offset = within_block
        return virtual_offset

    def read(self, size=-1):
        if size < 0:
            raise NotImplementedError("Don't be greedy, that could be massive!")
        elif size == 0:
            return ""
        elif self._within_block_offset + size < len(self._buffer):
            data = self._buffer[self._within_block_offset:self._within_block_offset + size]
            self._within_block_offset += size
            return data
        else:
            data = self._buffer[self._within_block_offset:]
            size -= len(data)
            self._load_block() #will reset offsets
            if not self._buffer:
                return data #EOF
            else:
                #TODO - Avoid recursion
                return data + self.read(size)

    def readline(self):
        i = self._buffer.find("\n", self._within_block_offset)
        if i != -1:
            data = self._buffer[self._within_block_offset:i+1]
            self._within_block_offset = i + 1
            assert data.endswith("\n")
            return data
        else:
            data = self._buffer[self._within_block_offset:]
            self._load_block() #will reset offsets
            if not self._buffer:
                return data #EOF
            else:
                #TODO - Avoid recursion
                return data + self.readline()

    def close(self):
        self._handle.close()
        self._buffer = None
        self._block_start_offset = None


class BgzfWriter(object):

    def __init__(self, filename=None, mode=None, fileobj=None, compresslevel=6):
        if fileobj:
            assert filename is None and mode is None            
            handle = fileobj
        else:
            if "w" not in mode.lower() \
            and "a" not in mode.lower():
                raise ValueError("Must use write or append mode, not %r" % mode)
            handle = __builtin__.open(filename, mode)
        self._handle = handle
        self._buffer = "" #Bytes string
        self.compresslevel = compresslevel

    def _write_block(self, block):
        #print "Saving %i bytes" % len(block)
        start_offset = self._handle.tell()
        assert len(block) <= 65536
        #Giving a negative window bits means no gzip/zlib headers, -15 used in samtools
        c = zlib.compressobj(self.compresslevel,
                             zlib.DEFLATED,
                             -15,
                             zlib.DEF_MEM_LEVEL,
                             0)
        compressed = c.compress(block) + c.flush()
        del c
        assert len(compressed) < 65536, "TODO - Didn't compress enough, try less data in this block"
        crc = zlib.crc32(block)
        #Should cope with a mix of Python platforms...
        if crc < 0:
            crc = struct.pack("<i", crc)
        else:
            crc = struct.pack("<I", crc)
        bsize = struct.pack("<H", len(compressed)+25) #includes -1
        crc = struct.pack("<I", zlib.crc32(block) & 0xffffffffL)
        uncompressed_length = struct.pack("<I", len(block))
        #Fixed 16 bytes,
        # gzip magic bytes (4) mod time (4),
        # gzip flag (1), os (1), extra length which is six (2),
        # sub field which is BC (2), sub field length of two (2),
        #Variable data,
        #2 bytes: block length as BC sub field (2)
        #X bytes: the data
        #8 bytes: crc (4), uncompressed data length (4)
        data = "\x1f\x8b\x08\x04\x00\x00\x00\x00" \
             + "\x00\xff\x06\x00\x42\x43\x02\x00" \
             + bsize + compressed + crc + uncompressed_length
        self._handle.write(data)

    def write(self, data):
        #block_size = 2**16 = 65536
        data_len = len(data)
        if len(self._buffer) + data_len < 65536:
            #print "Cached %r" % data
            self._buffer += data
            return
        else:
            #print "Got %r, writing out some data..." % data
            self._buffer += data
            while len(self._buffer) >= 65536:
                self._write_block(self._buffer[:65536])
                self._buffer = self._buffer[65536:]

    def flush(self):
        while len(self._buffer) >= 65536:
            self._write_block(self._buffer[:65535])
            self._buffer = self._buffer[65535:]
        self._write_block(self._buffer)
        self._buffer = ""
        self._handle.flush()

    def close(self):
        if self._buffer:
            self.flush()
        #samtools seems to look for a magic EOF marker, just a 28 byte empty BGZF block,
        #but we'll leave it up to the calling code to do that (e.g. flush twice).
        #self._handle.write("\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC")
        #self._handle.write("\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00")
        #self._handle.flush()
        self._handle.close()

    def tell(self):
        """Returns a BGZF 64-bit virtual offset."""
        return make_virtual_offset(self.handle.tell(), len(self._buffer)) 

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        print "Call this with no arguments and pipe uncompressed data in on stdin"
        print "and it will produce BGZF compressed data on stdout. e.g."
        print
        print "./bgzf.py < example.fastq > example.fastq.bgz"
        print
        print "The extension convention of *.bgz is to distinugish these from *.gz"
        print "used for standard gzipped files without the block structure of BGZF."
        print "You can use the standard gunzip command to decompress BGZF files,"
        print "if it complains about the extension try something like this:"
        print
        print "cat example.fastq.bgz | gunzip > example.fastq"
        print
        print "See also the tool bgzip that comes with samtools"
        sys.exit(0)

    sys.stderr.write("Producing BGZF output from stdin...\n")
    w = BgzfWriter(fileobj=sys.stdout)
    while True:
        data = sys.stdin.read(65536)
        w.write(data)
        if not data:
            break
    #Works but much slower
    #for data in sys.stdin:
    #    w.write(data)
    w.flush()
    w.flush() #Double flush triggers an empty BGZF block as EOF marker
    w.close()
    sys.stderr.write("BGZF data produced\n")
