# Copyright 2020 by Michiel de Hoon
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Bio.SeqIO support for UCSC's "twoBit" (.2bit) file format."""

# The .2bit file format is defined by UCSC as follows
# (see http://genome.ucsc.edu/FAQ/FAQformat.html#format7):
#
# ----------------------------------------------------------------------------
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
# ----------------------------------------------------------------------------
#
# This parser reads the index stored in the twoBit file, as well as the masked
# regions and the N's for ean sequence. It also creates sequence data objects
# (TwoBitSequence objects, defined in Bio/SeqIO/_twoBitIO.c), which support only
# two methods: __len__ and __getitem__. The former will return the length of the
# sequence, while the latter returns the sequence (as a bytes object) for the
# requested region.  Using the information in the index, The __getitem__ method
# calculates the file position at which the requested region starts, and only
# reads the requested sequence region. Note that the full sequence of a record
# is loaded only if specifically requested, making the parser memory-efficient.
#
# The TwoBitIterator object implements the __getitem__, keys, and __len__
# methods that allow it to be used as a dictionary.


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import _twoBitIO
from .Interfaces import SequenceIterator


class TwoBitIterator(SequenceIterator):
    """Parser for UCSC twoBit (.2bit) files."""

    def __init__(self, source):
        """Read the file index."""
        super().__init__(source, mode="b", fmt="twoBit")
        isByteSwapped, names, sequences = _twoBitIO.TwoBitIterator(self.stream)
        self.isByteSwapped = isByteSwapped
        self.names = [name.decode("ASCII") for name in names]
        self.sequences = [Seq(sequence) for sequence in sequences]
        # wait to close the file until the TwoBitIterator goes out of scope:
        self.should_close_stream = False

    def parse(self, stream):
        """Iterate over the sequences in the file."""
        for name, sequence in zip(self.names, self.sequences):
            sequence = Seq(sequence)
            record = SeqRecord(sequence, id=name)
            yield record

    def __getitem__(self, key):
        try:
            index = self.names.index(key)
        except ValueError:
            raise KeyError(key) from None
        sequence = self.sequences[index]
        return SeqRecord(sequence, id=key)

    def keys(self):
        """Return a list with the names of the sequences in the file."""
        return self.names

    def __len__(self):
        return len(self.names)
