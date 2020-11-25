# Copyright 2020 by Michiel de Hoon
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Bio.SeqIO support for UCSC's "twoBit" (.2bit) file format.

This parser reads the index stored in the twoBit file, as well as the masked
regions and the N's for ean sequence. It also creates sequence data objects
(TwoBitSequenceData objects), which support only two methods: __len__ and
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
from Bio.SeqRecord import SeqRecord

from .Interfaces import SequenceIterator


class TwoBitSequenceData():

    bases = numpy.array([b"TTTT",  # 00 00 00 00
                         b"TTTC",  # 00 00 00 01
                         b"TTTA",  # 00 00 00 10
                         b"TTTG",  # 00 00 00 11
                         b"TTCT",  # 00 00 01 00
                         b"TTCC",  # 00 00 01 01
                         b"TTCA",  # 00 00 01 10
                         b"TTCG",  # 00 00 01 11
                         b"TTAT",  # 00 00 10 00
                         b"TTAC",  # 00 00 10 01
                         b"TTAA",  # 00 00 10 10
                         b"TTAG",  # 00 00 10 11
                         b"TTGT",  # 00 00 11 00
                         b"TTGC",  # 00 00 11 01
                         b"TTGA",  # 00 00 11 10
                         b"TTGG",  # 00 00 11 11
                         b"TCTT",  # 00 01 00 00
                         b"TCTC",  # 00 01 00 01
                         b"TCTA",  # 00 01 00 10
                         b"TCTG",  # 00 01 00 11
                         b"TCCT",  # 00 01 01 00
                         b"TCCC",  # 00 01 01 01
                         b"TCCA",  # 00 01 01 10
                         b"TCCG",  # 00 01 01 11
                         b"TCAT",  # 00 01 10 00
                         b"TCAC",  # 00 01 10 01
                         b"TCAA",  # 00 01 10 10
                         b"TCAG",  # 00 01 10 11
                         b"TCGT",  # 00 01 11 00
                         b"TCGC",  # 00 01 11 01
                         b"TCGA",  # 00 01 11 10
                         b"TCGG",  # 00 01 11 11
                         b"TATT",  # 00 10 00 00
                         b"TATC",  # 00 10 00 01
                         b"TATA",  # 00 10 00 10
                         b"TATG",  # 00 10 00 11
                         b"TACT",  # 00 10 01 00
                         b"TACC",  # 00 10 01 01
                         b"TACA",  # 00 10 01 10
                         b"TACG",  # 00 10 01 11
                         b"TAAT",  # 00 10 10 00
                         b"TAAC",  # 00 10 10 01
                         b"TAAA",  # 00 10 10 10
                         b"TAAG",  # 00 10 10 11
                         b"TAGT",  # 00 10 11 00
                         b"TAGC",  # 00 10 11 01
                         b"TAGA",  # 00 10 11 10
                         b"TAGG",  # 00 10 11 11
                         b"TGTT",  # 00 11 00 00
                         b"TGTC",  # 00 11 00 01
                         b"TGTA",  # 00 11 00 10
                         b"TGTG",  # 00 11 00 11
                         b"TGCT",  # 00 11 01 00
                         b"TGCC",  # 00 11 01 01
                         b"TGCA",  # 00 11 01 10
                         b"TGCG",  # 00 11 01 11
                         b"TGAT",  # 00 11 10 00
                         b"TGAC",  # 00 11 10 01
                         b"TGAA",  # 00 11 10 10
                         b"TGAG",  # 00 11 10 11
                         b"TGGT",  # 00 11 11 00
                         b"TGGC",  # 00 11 11 01
                         b"TGGA",  # 00 11 11 10
                         b"TGGG",  # 00 11 11 11
                         b"CTTT",  # 01 00 00 00
                         b"CTTC",  # 01 00 00 01
                         b"CTTA",  # 01 00 00 10
                         b"CTTG",  # 01 00 00 11
                         b"CTCT",  # 01 00 01 00
                         b"CTCC",  # 01 00 01 01
                         b"CTCA",  # 01 00 01 10
                         b"CTCG",  # 01 00 01 11
                         b"CTAT",  # 01 00 10 00
                         b"CTAC",  # 01 00 10 01
                         b"CTAA",  # 01 00 10 10
                         b"CTAG",  # 01 00 10 11
                         b"CTGT",  # 01 00 11 00
                         b"CTGC",  # 01 00 11 01
                         b"CTGA",  # 01 00 11 10
                         b"CTGG",  # 01 00 11 11
                         b"CCTT",  # 01 01 00 00
                         b"CCTC",  # 01 01 00 01
                         b"CCTA",  # 01 01 00 10
                         b"CCTG",  # 01 01 00 11
                         b"CCCT",  # 01 01 01 00
                         b"CCCC",  # 01 01 01 01
                         b"CCCA",  # 01 01 01 10
                         b"CCCG",  # 01 01 01 11
                         b"CCAT",  # 01 01 10 00
                         b"CCAC",  # 01 01 10 01
                         b"CCAA",  # 01 01 10 10
                         b"CCAG",  # 01 01 10 11
                         b"CCGT",  # 01 01 11 00
                         b"CCGC",  # 01 01 11 01
                         b"CCGA",  # 01 01 11 10
                         b"CCGG",  # 01 01 11 11
                         b"CATT",  # 01 10 00 00
                         b"CATC",  # 01 10 00 01
                         b"CATA",  # 01 10 00 10
                         b"CATG",  # 01 10 00 11
                         b"CACT",  # 01 10 01 00
                         b"CACC",  # 01 10 01 01
                         b"CACA",  # 01 10 01 10
                         b"CACG",  # 01 10 01 11
                         b"CAAT",  # 01 10 10 00
                         b"CAAC",  # 01 10 10 01
                         b"CAAA",  # 01 10 10 10
                         b"CAAG",  # 01 10 10 11
                         b"CAGT",  # 01 10 11 00
                         b"CAGC",  # 01 10 11 01
                         b"CAGA",  # 01 10 11 10
                         b"CAGG",  # 01 10 11 11
                         b"CGTT",  # 01 11 00 00
                         b"CGTC",  # 01 11 00 01
                         b"CGTA",  # 01 11 00 10
                         b"CGTG",  # 01 11 00 11
                         b"CGCT",  # 01 11 01 00
                         b"CGCC",  # 01 11 01 01
                         b"CGCA",  # 01 11 01 10
                         b"CGCG",  # 01 11 01 11
                         b"CGAT",  # 01 11 10 00
                         b"CGAC",  # 01 11 10 01
                         b"CGAA",  # 01 11 10 10
                         b"CGAG",  # 01 11 10 11
                         b"CGGT",  # 01 11 11 00
                         b"CGGC",  # 01 11 11 01
                         b"CGGA",  # 01 11 11 10
                         b"CGGG",  # 01 11 11 11
                         b"ATTT",  # 10 00 00 00
                         b"ATTC",  # 10 00 00 01
                         b"ATTA",  # 10 00 00 10
                         b"ATTG",  # 10 00 00 11
                         b"ATCT",  # 10 00 01 00
                         b"ATCC",  # 10 00 01 01
                         b"ATCA",  # 10 00 01 10
                         b"ATCG",  # 10 00 01 11
                         b"ATAT",  # 10 00 10 00
                         b"ATAC",  # 10 00 10 01
                         b"ATAA",  # 10 00 10 10
                         b"ATAG",  # 10 00 10 11
                         b"ATGT",  # 10 00 11 00
                         b"ATGC",  # 10 00 11 01
                         b"ATGA",  # 10 00 11 10
                         b"ATGG",  # 10 00 11 11
                         b"ACTT",  # 10 01 00 00
                         b"ACTC",  # 10 01 00 01
                         b"ACTA",  # 10 01 00 10
                         b"ACTG",  # 10 01 00 11
                         b"ACCT",  # 10 01 01 00
                         b"ACCC",  # 10 01 01 01
                         b"ACCA",  # 10 01 01 10
                         b"ACCG",  # 10 01 01 11
                         b"ACAT",  # 10 01 10 00
                         b"ACAC",  # 10 01 10 01
                         b"ACAA",  # 10 01 10 10
                         b"ACAG",  # 10 01 10 11
                         b"ACGT",  # 10 01 11 00
                         b"ACGC",  # 10 01 11 01
                         b"ACGA",  # 10 01 11 10
                         b"ACGG",  # 10 01 11 11
                         b"AATT",  # 10 10 00 00
                         b"AATC",  # 10 10 00 01
                         b"AATA",  # 10 10 00 10
                         b"AATG",  # 10 10 00 11
                         b"AACT",  # 10 10 01 00
                         b"AACC",  # 10 10 01 01
                         b"AACA",  # 10 10 01 10
                         b"AACG",  # 10 10 01 11
                         b"AAAT",  # 10 10 10 00
                         b"AAAC",  # 10 10 10 01
                         b"AAAA",  # 10 10 10 10
                         b"AAAG",  # 10 10 10 11
                         b"AAGT",  # 10 10 11 00
                         b"AAGC",  # 10 10 11 01
                         b"AAGA",  # 10 10 11 10
                         b"AAGG",  # 10 10 11 11
                         b"AGTT",  # 10 11 00 00
                         b"AGTC",  # 10 11 00 01
                         b"AGTA",  # 10 11 00 10
                         b"AGTG",  # 10 11 00 11
                         b"AGCT",  # 10 11 01 00
                         b"AGCC",  # 10 11 01 01
                         b"AGCA",  # 10 11 01 10
                         b"AGCG",  # 10 11 01 11
                         b"AGAT",  # 10 11 10 00
                         b"AGAC",  # 10 11 10 01
                         b"AGAA",  # 10 11 10 10
                         b"AGAG",  # 10 11 10 11
                         b"AGGT",  # 10 11 11 00
                         b"AGGC",  # 10 11 11 01
                         b"AGGA",  # 10 11 11 10
                         b"AGGG",  # 10 11 11 11
                         b"GTTT",  # 11 00 00 00
                         b"GTTC",  # 11 00 00 01
                         b"GTTA",  # 11 00 00 10
                         b"GTTG",  # 11 00 00 11
                         b"GTCT",  # 11 00 01 00
                         b"GTCC",  # 11 00 01 01
                         b"GTCA",  # 11 00 01 10
                         b"GTCG",  # 11 00 01 11
                         b"GTAT",  # 11 00 10 00
                         b"GTAC",  # 11 00 10 01
                         b"GTAA",  # 11 00 10 10
                         b"GTAG",  # 11 00 10 11
                         b"GTGT",  # 11 00 11 00
                         b"GTGC",  # 11 00 11 01
                         b"GTGA",  # 11 00 11 10
                         b"GTGG",  # 11 00 11 11
                         b"GCTT",  # 11 01 00 00
                         b"GCTC",  # 11 01 00 01
                         b"GCTA",  # 11 01 00 10
                         b"GCTG",  # 11 01 00 11
                         b"GCCT",  # 11 01 01 00
                         b"GCCC",  # 11 01 01 01
                         b"GCCA",  # 11 01 01 10
                         b"GCCG",  # 11 01 01 11
                         b"GCAT",  # 11 01 10 00
                         b"GCAC",  # 11 01 10 01
                         b"GCAA",  # 11 01 10 10
                         b"GCAG",  # 11 01 10 11
                         b"GCGT",  # 11 01 11 00
                         b"GCGC",  # 11 01 11 01
                         b"GCGA",  # 11 01 11 10
                         b"GCGG",  # 11 01 11 11
                         b"GATT",  # 11 10 00 00
                         b"GATC",  # 11 10 00 01
                         b"GATA",  # 11 10 00 10
                         b"GATG",  # 11 10 00 11
                         b"GACT",  # 11 10 01 00
                         b"GACC",  # 11 10 01 01
                         b"GACA",  # 11 10 01 10
                         b"GACG",  # 11 10 01 11
                         b"GAAT",  # 11 10 10 00
                         b"GAAC",  # 11 10 10 01
                         b"GAAA",  # 11 10 10 10
                         b"GAAG",  # 11 10 10 11
                         b"GAGT",  # 11 10 11 00
                         b"GAGC",  # 11 10 11 01
                         b"GAGA",  # 11 10 11 10
                         b"GAGG",  # 11 10 11 11
                         b"GGTT",  # 11 11 00 00
                         b"GGTC",  # 11 11 00 01
                         b"GGTA",  # 11 11 00 10
                         b"GGTG",  # 11 11 00 11
                         b"GGCT",  # 11 11 01 00
                         b"GGCC",  # 11 11 01 01
                         b"GGCA",  # 11 11 01 10
                         b"GGCG",  # 11 11 01 11
                         b"GGAT",  # 11 11 10 00
                         b"GGAC",  # 11 11 10 01
                         b"GGAA",  # 11 11 10 10
                         b"GGAG",  # 11 11 10 11
                         b"GGGT",  # 11 11 11 00
                         b"GGGC",  # 11 11 11 01
                         b"GGGA",  # 11 11 11 10
                         b"GGGG",  # 11 11 11 11
                        ], dtype=bytes)
    bases.dtype = "uint8"
    bases.shape = (256, 4)

    def __init__(self, stream, offset):
        self.stream = stream
        self.offset = offset

    def __getitem__(self, key):
        length = self.length
        if isinstance(key, slice):
            start, end, step = key.indices(length)
            size = len(range(start, end, step))
            if size == 0:
                return b""
        else:
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
        sequence = self.bases[data]
        sequence.shape = (4*byteSize, )
        skip = byteStart * 4
        sequence = sequence[start-skip:start-skip+size]
        nBlockStarts = self.nBlockStarts
        nBlockEnds = self.nBlockEnds
        i = nBlockEnds.searchsorted(start+1)
        j = i + nBlockStarts[i:].searchsorted(end)
        n = ord("N")
        for k in range(i, j):
            nBlockStart = nBlockStarts[k]
            nBlockEnd = nBlockEnds[k]
            if nBlockStart < start:
                nBlockStart = start
            if end < nBlockEnd:
                nBlockEnd = end
            sequence[nBlockStart - start : nBlockEnd - start] = n
        maskBlockStarts = self.maskBlockStarts
        maskBlockEnds = self.maskBlockEnds
        i = maskBlockEnds.searchsorted(start+1)
        j = i + maskBlockStarts[i:].searchsorted(end)
        difference = ord("a") - ord("A")
        for k in range(i, j):
            maskBlockStart = maskBlockStarts[k]
            maskBlockEnd = maskBlockEnds[k]
            if maskBlockStart < start:
                maskBlockStart = start
            if end < maskBlockEnd:
                maskBlockEnd = end
            sequence[maskBlockStart - start : maskBlockEnd - start] += difference
        if step == 1:
            return bytes(sequence)
        else:
            return bytes(sequence[::step])

    def __len__(self):
        return self.length


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
                "version-1 twoBit files with 64-bit offsets for index are currently not supported")
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
            sequence = TwoBitSequenceData(stream, offset)
            sequences[name] = sequence
        self.sequences = sequences
        for name, sequence in sequences.items():
            offset = sequence.offset
            stream.seek(offset)
            data = stream.read(4)
            dnaSize = int.from_bytes(data, byteorder, signed=False)
            sequence.length = dnaSize
            data = stream.read(4)
            nBlockCount = int.from_bytes(data, byteorder, signed=False)
            nBlockStarts = numpy.fromfile(stream, dtype=dtype, count=nBlockCount)
            nBlockSizes = numpy.fromfile(stream, dtype=dtype, count=nBlockCount)
            sequence.nBlockStarts = nBlockStarts
            sequence.nBlockEnds = nBlockStarts + nBlockSizes
            data = stream.read(4)
            maskBlockCount = int.from_bytes(data, byteorder, signed=False)
            maskBlockStarts = numpy.fromfile(stream, dtype=dtype, count=maskBlockCount)
            maskBlockSizes = numpy.fromfile(stream, dtype=dtype, count=maskBlockCount)
            sequence.maskBlockStarts = maskBlockStarts
            sequence.maskBlockEnds = maskBlockStarts + maskBlockSizes
            data = stream.read(4)
            reserved = int.from_bytes(data, byteorder, signed=False)
            if reserved != 0:
                raise ValueError("Found non-zero reserved field %u" % reserved)
            sequence.offset = stream.tell()

    def parse(self, stream):
        """Iterate over the sequences in the file."""
        for name, sequence in self.sequences.items():
            sequence = Seq(sequence)
            record = SeqRecord(sequence, id=name)
            yield record

    def __getitem__(self, name):
        try:
            sequence = self.sequences[name]
        except ValueError:
            raise KeyError(key) from None
        sequence = Seq(sequence)
        return SeqRecord(sequence, id=name)

    def keys(self):
        """Return a list with the names of the sequences in the file."""
        return self.sequences.keys()

    def __len__(self):
        return len(self.sequences)
