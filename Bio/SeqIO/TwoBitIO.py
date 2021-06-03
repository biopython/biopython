# Copyright 2020 by Michiel de Hoon
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SeqIO support for UCSC's "twoBit" (.2bit) file format.

This parser reads the index stored in the twoBit file, as well as the masked
regions and the N's for each sequence. It also creates sequence data objects
(_TwoBitSequenceData objects), which support only two methods: __len__ and
__getitem__. The former will return the length of the sequence, while the
latter returns the sequence (as a bytes object) for the requested region.

Using the information in the index, the __getitem__ method calculates the file
position at which the requested region starts, and only reads the requested
sequence region. Note that the full sequence of a record is loaded only if
specifically requested, making the parser memory-efficient.

The TwoBitIterator object implements the __getitem__, keys, and __len__
methods that allow it to be used as a dictionary.
"""
# The .2bit file format is defined by UCSC as follows
# (see http://genome.ucsc.edu/FAQ/FAQformat.html#format7):
#
#
# A .2bit file stores multiple DNA sequences (up to 4 Gb total) in a compact
# randomly-accessible format. The file contains masking information as well
# as the DNA itself.
#
# The file begins with a 16-byte header containing the following fields:
#
# signature - the number 0x1A412743 in the architecture of the machine that
#             created the file
# version - zero for now. Readers should abort if they see a version number
#           higher than 0
# sequenceCount - the number of sequences in the file
# reserved - always zero for now
#
# All fields are 32 bits unless noted. If the signature value is not as
# given, the reader program should byte-swap the signature and check if the
# swapped version matches. If so, all multiple-byte entities in the file
# will have to be byte-swapped. This enables these binary files to be used
# unchanged on different architectures.
#
# The header is followed by a file index, which contains one entry for each
# sequence. Each index entry contains three fields:
#
# nameSize - a byte containing the length of the name field
# name - the sequence name itself (in ASCII-compatible byte string), of
#        variable length depending on nameSize
# offset - the 32-bit offset of the sequence data relative to the start of
#          the file, not aligned to any 4-byte padding boundary
#
# The index is followed by the sequence records, which contain nine fields:
#
# dnaSize - number of bases of DNA in the sequence
# nBlockCount - the number of blocks of Ns in the file (representing unknown
#               sequence)
# nBlockStarts - an array of length nBlockCount of 32 bit integers
#                indicating the (0-based) starting position of a block of Ns
# nBlockSizes - an array of length nBlockCount of 32 bit integers indicating
#               the length of a block of Ns
# maskBlockCount - the number of masked (lower-case) blocks
# maskBlockStarts - an array of length maskBlockCount of 32 bit integers
#                   indicating the (0-based) starting position of a masked block
# maskBlockSizes - an array of length maskBlockCount of 32 bit integers
#                  indicating the length of a masked block
# reserved - always zero for now
# packedDna - the DNA packed to two bits per base, represented as so:
#             T - 00, C - 01, A - 10, G - 11. The first base is in the most
#             significant 2-bit byte; the last base is in the least significan
#             2 bits. For example, the sequence TCAG is represented as 00011011.
import numpy

from Bio.Seq import Seq
from Bio.Seq import SequenceDataAbstractBaseClass
from Bio.SeqRecord import SeqRecord

from . import _twoBitIO
from .Interfaces import SequenceIterator


class _TwoBitSequenceData(SequenceDataAbstractBaseClass):
    """Stores information needed to retrieve sequence data from a .2bit file (PRIVATE).

    Objects of this class store the file position at which the sequence data
    start, the sequence length, and the start and end position of unknown (N)
    and masked (lowercase) letters in the sequence.

    Only two methods are provided: __len__ and __getitem__. The former will
    return the length of the sequence, while the latter returns the sequence
    (as a bytes object) for the requested region. The full sequence of a record
    is loaded only if explicitly requested.
    """

    __slots__ = ("stream", "offset", "length", "nBlocks", "maskBlocks")

    def __init__(self, stream, offset, length):
        """Initialize the file stream and file position of the sequence data."""
        self.stream = stream
        self.offset = offset
        self.length = length
        super().__init__()

    def __getitem__(self, key):
        """Return the sequence contents (as a bytes object) for the requested region."""
        length = self.length
        if isinstance(key, slice):
            start, end, step = key.indices(length)
            size = len(range(start, end, step))
            if size == 0:
                return b""
        else:
            if key < 0:
                key += length
                if key < 0:
                    raise IndexError("index out of range")
            start = key
            end = key + 1
            step = 1
            size = 1
        byteStart = start // 4
        byteEnd = (end + 3) // 4
        byteSize = byteEnd - byteStart
        stream = self.stream
        try:
            stream.seek(self.offset + byteStart)
        except ValueError as exception:
            if str(exception) == "seek of closed file":
                raise ValueError("cannot retrieve sequence: file is closed") from None
            raise
        data = numpy.fromfile(stream, dtype="uint8", count=byteSize)
        sequence = _twoBitIO.convert(
            data, start, end, step, self.nBlocks, self.maskBlocks
        )
        if isinstance(key, slice):
            return sequence
        else:  # single nucleotide
            return ord(sequence)

    def __len__(self):
        """Get the sequence length."""
        return self.length

    def upper(self):
        """Remove the sequence mask."""
        data = _TwoBitSequenceData(self.stream, self.offset, self.length)
        data.nBlocks = self.nBlocks[:, :]
        data.maskBlocks = numpy.empty((0, 2), dtype="uint32")
        return data

    def lower(self):
        """Extend the sequence mask to the full sequence."""
        data = _TwoBitSequenceData(self.stream, self.offset, self.length)
        data.nBlocks = self.nBlocks[:, :]
        data.maskBlocks = numpy.array([[0, self.length]], dtype="uint32")
        return data


class TwoBitIterator(SequenceIterator):
    """Parser for UCSC twoBit (.2bit) files."""

    def __init__(self, source):
        """Read the file index."""
        super().__init__(source, mode="b", fmt="twoBit")
        # wait to close the file until the TwoBitIterator goes out of scope:
        self.should_close_stream = False
        stream = self.stream
        data = stream.read(4)
        if not data:
            raise ValueError("Empty file.")
        byteorders = ("little", "big")
        dtypes = ("<u4", ">u4")
        for byteorder, dtype in zip(byteorders, dtypes):
            signature = int.from_bytes(data, byteorder)
            if signature == 0x1A412743:
                break
        else:
            raise ValueError("Unknown signature")
        self.byteorder = byteorder
        data = stream.read(4)
        version = int.from_bytes(data, byteorder, signed=False)
        if version == 1:
            raise ValueError(
                "version-1 twoBit files with 64-bit offsets for index are currently not supported"
            )
        if version != 0:
            raise ValueError("Found unexpected file version %u; aborting" % version)
        data = stream.read(4)
        sequenceCount = int.from_bytes(data, byteorder, signed=False)
        data = stream.read(4)
        reserved = int.from_bytes(data, byteorder, signed=False)
        if reserved != 0:
            raise ValueError("Found non-zero reserved field; aborting")
        sequences = {}
        for i in range(sequenceCount):
            data = stream.read(1)
            nameSize = int.from_bytes(data, byteorder, signed=False)
            data = stream.read(nameSize)
            name = data.decode("ASCII")
            data = stream.read(4)
            offset = int.from_bytes(data, byteorder, signed=False)
            sequences[name] = (stream, offset)
        self.sequences = sequences
        for name, (stream, offset) in sequences.items():
            stream.seek(offset)
            data = stream.read(4)
            dnaSize = int.from_bytes(data, byteorder, signed=False)
            sequence = _TwoBitSequenceData(stream, offset, dnaSize)
            data = stream.read(4)
            nBlockCount = int.from_bytes(data, byteorder, signed=False)
            nBlockStarts = numpy.fromfile(stream, dtype=dtype, count=nBlockCount)
            nBlockSizes = numpy.fromfile(stream, dtype=dtype, count=nBlockCount)
            sequence.nBlocks = numpy.empty((nBlockCount, 2), dtype="uint32")
            sequence.nBlocks[:, 0] = nBlockStarts
            sequence.nBlocks[:, 1] = nBlockStarts + nBlockSizes
            data = stream.read(4)
            maskBlockCount = int.from_bytes(data, byteorder, signed=False)
            maskBlockStarts = numpy.fromfile(stream, dtype=dtype, count=maskBlockCount)
            maskBlockSizes = numpy.fromfile(stream, dtype=dtype, count=maskBlockCount)
            sequence.maskBlocks = numpy.empty((maskBlockCount, 2), dtype="uint32")
            sequence.maskBlocks[:, 0] = maskBlockStarts
            sequence.maskBlocks[:, 1] = maskBlockStarts + maskBlockSizes
            data = stream.read(4)
            reserved = int.from_bytes(data, byteorder, signed=False)
            if reserved != 0:
                raise ValueError("Found non-zero reserved field %u" % reserved)
            sequence.offset = stream.tell()
            sequences[name] = sequence

    def parse(self, stream):
        """Iterate over the sequences in the file."""
        for name, sequence in self.sequences.items():
            sequence = Seq(sequence)
            record = SeqRecord(sequence, id=name)
            yield record

    def __getitem__(self, name):
        """Return sequence associated with given name as a SeqRecord object."""
        try:
            sequence = self.sequences[name]
        except ValueError:
            raise KeyError(name) from None
        sequence = Seq(sequence)
        return SeqRecord(sequence, id=name)

    def keys(self):
        """Return a list with the names of the sequences in the file."""
        return self.sequences.keys()

    def __len__(self):
        """Return number of sequences."""
        return len(self.sequences)
