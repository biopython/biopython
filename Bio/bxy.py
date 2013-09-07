# Copyright 2012-2013 by Peter Cock.
# All rights reserved.
"""Python module for random access to blocked XY files (LMZA compression).

Python 3.3 includes the moudle lmza which allows read and write access
to XY files, however random access is emulated and is slow (much like
with Python's gzip module).

Biopython's bgzf module supports efficient random access to BGZF files,
Blocked GNU Zip Format, a variant of GZIP using blocks for random access.
This was implemented as a Python module calling the library zlib (in C)
for efficient compression and decompression. Here we take a similar
approach.
"""

import os
import lzma
import struct
from io import BytesIO
from zlib import crc32

from Bio._py3k import _as_bytes, _as_string

_stream_header_magic = _as_bytes("\xfd7zXZ\x00")
_stream_footer_magic = _as_bytes("YZ")
_empty_bytes_string = _as_bytes("")
_null = _as_bytes("\0")

def _get_variable_int(handle):
    """Read a variable length encoded integer from XZ file handle.

    Returns a tuple, number of bytes used, and the value.
    
    Based on the C function decode given in the XY specification.
    Essentially this uses base 128 (i.e. 2**7), with the top bit
    set to indicate another byte/digit is required, giving the
    least significant digit/byte first.

    TODO - Make this work from a bytes string & offset
    (can return the new offset and the value).
    """
    #TODO - Is this endian safe?
    i = 0
    b = handle.read(1)[0] #As an integer
    value = b & 0x7F #Mask out the high bit (if set)
    buf = [b]
    #print(i, b, value)
    while b & 0x80: #If high bit 0x80 = 128 is set
        i += 1
        b = handle.read(1)[0] #As an integer
        buf.append(b)
        if (i >= 9) or (b == _null):
            raise OverFlowError
        #Add new value is (b & 0x7F) multiplied by 128**i
        #which we can calculate using an i*7 shift
        #value |= (b & 0x7F) << (i * 7)
        value += (b & 0x7F)*(128**i)
        #print(i, b, value)
    return i+1, value

#Quick self test:
assert (1, 7) == _get_variable_int(BytesIO(_as_bytes('\x07')))
assert (3, 1024**2) == _get_variable_int(BytesIO(_as_bytes('\x80\x80@')))

def _load_stream_header(handle):
    header = handle.read(12)
    assert header[0:6] == _stream_header_magic, header[0:6]
    stream_flags, crc = struct.unpack("<HI", header[6:12])
    assert crc32(header[6:8]) == crc #Checksum for stream flags
    return stream_flags, crc

#_footer_fmt =
def _load_stream_footer(handle):
    footer = handle.read(12)
    assert footer[10:12] == _stream_footer_magic, footer[10:12]
    crc, stored_backward_size, stream_flags = struct.unpack("<IIH", footer[0:10])
    #Checksum for backward size and stream flags:
    assert crc32(footer[4:10]) == crc
    real_backward_size = (stored_backward_size + 1) * 4;
    #TODO - Check for null padding?
    return crc, real_backward_size, stream_flags

def _load_stream_block_header(handle, expected_comp_size, expected_uncomp_size):
    #print("Loading block header, expect comp %i --> %i" % (expected_uncomp_size, expected_comp_size))
    block_header_size, block_flags = struct.unpack("<BB", handle.read(2))
    #Can't be zero! 
    assert 0x01 < block_header_size <= 0xFF, block_header_size
    filter_count = (block_flags & 0x03) + 1 #allowed 1 to 4 filters
    if block_flags & 0x3C:
        raise ValueError("Bits 0x3c should be zero in block flags, got %i (0x%x)" % (block_flags, block_flags))
    if block_flags & 0x40:
        bytes_used, compressed_size = _get_variable_int(handle)
        assert compressed_size == expected_comp_size, \
            "Compressed size in block header %i, expected %i" \
            % (compressed_size, expected_comp_size)
    if block_flags & 0x80:
        bytes_used, uncompressed_size = _get_variable_int(handle)
        assert uncompressed_size == expected_uncomp_size, \
            "Uncompressed size in block header %i, expected %i" \
            % (uncompressed_size, expected_uncomp_size)
    #TODO - Filters, padding, CRC32

def _load_stream_index(handle, expected_size):
    indicator = handle.read(1)
    assert indicator == _null, indicator
    bytes_used, record_count = _get_variable_int(handle)
    size = 1 + bytes_used
    #print("Index contains %i records" % record_count)
    for record in range(record_count):
        bytes_used, unpadded_size = _get_variable_int(handle)
        size += bytes_used
        bytes_used, uncompressed_size = _get_variable_int(handle)
        size += bytes_used
        #print("Record %i, unpadded size %i, uncompressed size %i" % (record, unpadded_size, uncompressed_size))
        yield unpadded_size, uncompressed_size
    if size % 4 != 0:
        #Null padding...
        bytes_used = 4 - (size % 4)
        padding = handle.read(bytes_used)
        assert padding == _null * bytes_used, padding
        size += bytes_used
    crc, = struct.unpack("<I", handle.read(4))
    #TODO - Checksum calculated for everything except the CRC
    size += 4
    assert size == expected_size, "Stream index size %i, expected %i" % (size, expected_size)


def _load_index(h):
    h.seek(0)
    _load_stream_header(h)
    size = os.fstat(h.fileno()).st_size
    print("File size is %i" % size)

    stream_end = size
    while stream_end > 0:
        print("---")
        print("Scanning stream ending at %i" % stream_end)
        h.seek(stream_end - 12)
        crc32, back_size, stream_flags = _load_stream_footer(h)
        h.seek(stream_end - 12 - back_size)
        block_sizes = list(_load_stream_index(h, back_size))
        total_comp_size = 0
        total_uncomp_size = 0
        for unpadded_size, uncompressed_size in block_sizes:
            #Round up unpadded size to multiple of four
            padded_size = unpadded_size + 4 - (unpadded_size % 4)
            #print("Block size %i (padded %i) --> %i" % (unpadded_size, padded_size, uncompressed_size))
            total_comp_size += padded_size
            total_uncomp_size += uncompressed_size
            #Sanity test, does the block header agree?
            block_start = stream_end - 12 - back_size - total_comp_size
            block_end = block_start + padded_size
            print("Block location %i to %i, size %i (padded %i) --> %i" \
                      % (block_start, block_end, unpadded_size, padded_size, uncompressed_size))
            h.seek(block_start)
            _load_stream_block_header(h, unpadded_size, uncompressed_size)
        print("Stream %i --> %i (%0.1f%%)" % (total_uncomp_size, total_comp_size, total_comp_size*100.0/total_uncomp_size))
        stream_start = size - 12 - back_size - total_comp_size - 12
        assert stream_start >= 0, stream_start
        h.seek(stream_start)
        print(_load_stream_header(h))
        print("Stream location %i to %i" % (stream_start, stream_end))
        #Check preceeding stream (if any)
        stream_end = stream_start
    print("--")
    assert stream_end == 0, stream_end


for f in ["thousand_blastx_nr.xml.xz",
          "thousand_blastx_nr.xml.b1M.xz",
          "thousand_blastx_nr.xml.s1M.xz"]:
    print("="*78)
    print(f)
    h = open(f, "rb")
    _load_index(h)
    print("Done")
    print()
