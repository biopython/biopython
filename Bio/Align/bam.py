"""Bio.Align support for the "bam" pairwise alignment format."""

import numpy as np
from Bio import bgzf
from io import StringIO, BytesIO
from typing import NamedTuple, Optional, Literal, Iterator
import struct
from Bio.Align import sam


def _to_uint(x):
    """Converts raw bytes in a BAM file to an unsigned integer"""
    return int.from_bytes(x, byteorder="little", signed=False)


class _ReferenceMetadata(NamedTuple):
    name: str
    length: int  # not used?


class _AlignmentInfo(NamedTuple):
    ref_id: int
    pos: int
    l_read_name: int
    mapq: int
    bin: int
    n_cigar_op: int
    flag: int
    l_seq: int
    next_ref_id: int
    next_pos: int
    tlen: int


class _BAMToSamReader(bgzf.BgzfReader):
    """Limited translation reader from BAM to SAM which allows seek(0), read(0), and iteration"""

    _references: Optional[list[_ReferenceMetadata]]
    _header_iterator: Optional[Iterator[str]]

    seq_encoding: bytes = bytes.maketrans(bytes(range(16)), b"=ACMGRSVTWYHKDBN")
    cigar_encoding: bytes = bytes.maketrans(b"012345678", b"MIDNSHP=X")
    qual_mapping: bytes = bytes(
        (letter + 33 if 0 <= letter < 94 else 255) for letter in range(256)
    )

    def __init__(self, filename=None, fileobj=None, max_cache=100):
        super().__init__(
            filename=filename, mode="rb", fileobj=fileobj, max_cache=max_cache
        )
        self._references = None
        self._header_iterator = None

    def seek(self, virtual_offset: Literal[0]) -> int:
        if virtual_offset != 0:
            raise NotImplementedError(
                f"BAM reader only accepts virtual_offset = 0, but found {virtual_offset=}"
            )
        self._header_iterator = None
        self._references = None
        return super().seek(virtual_offset)

    def read(self, size=-1):
        if size != 0:
            raise NotImplementedError(
                f"BAM reader only accepts size = 0, but found {size=}"
            )
        return ""

    @classmethod
    def seq_nibble_decode(cls, seq: bytes) -> bytes:
        seq_array = np.frombuffer(seq, dtype="<B")

        # split into high and low nibbles (this could probably be done faster at a low level if necessary)
        low, high = seq_array & 0x0F, seq_array >> 4

        # Interweave the high and low nibbles
        interweaved = np.ravel(np.column_stack((low, high)))
        return interweaved.tobytes().translate(cls.seq_encoding)

    @classmethod
    def cigar_nibble_decode(cls, cigar_bytes):
        cigar_arr = np.frombuffer(cigar_bytes, dtype="<I")

        # Split into high and low nibbles
        op, op_len = cigar_arr & 0x0F, cigar_arr >> 4

        # Convert to SAM operation characters
        translated_op = np.strings.translate(op.astype(np.bytes_), cls.cigar_encoding)
        # interweave op_len and op
        cigar = (np.char.array(op_len) + translated_op).tobytes().replace(b"\00", b"")
        return cigar

    def __next__(self):
        if self._header_iterator is None:
            magic_bam = super().read(4)
            assert magic_bam == b"BAM\x01", f"{magic_bam} != b'BAM\\x01'"
            header_length = _to_uint(super().read(4))
            self._header_iterator = iter(
                StringIO(super().read(header_length).decode("UTF-8"))
            )

        if self._references is None:
            # ----------------------------Header-------------------------------------
            try:
                return next(self._header_iterator)
            except StopIteration:
                pass

            # ----------------------------References-------------------------------------
            num_refs = _to_uint(super().read(4))
            self._references = []
            for i in range(num_refs):
                ref_name_length = _to_uint(super().read(4))

                ref_name = (
                    super().read(ref_name_length).rstrip(b"\00").decode("ASCII")
                )  # This should be ascii?
                ref_length = _to_uint(super().read(4))
                self._references.append(_ReferenceMetadata(ref_name, ref_length))

        # ----------------------------Alignments-------------------------------------
        block_size_bytes = super().read(4)
        if not block_size_bytes:
            raise StopIteration
        block_size = _to_uint(block_size_bytes)

        block = BytesIO(super().read(block_size))
        alignment_info = _AlignmentInfo._make(
            struct.unpack("<iiBBHHHIiii", block.read(32))
        )
        read_name = block.read(alignment_info.l_read_name).rstrip(b"\00")
        # TODO: handle case when there are too many operations
        cigar = (
            self.cigar_nibble_decode(block.read(alignment_info.n_cigar_op * 4))
            if alignment_info.n_cigar_op > 0
            else b"*"
        )

        seq = (
            self.seq_nibble_decode(block.read((alignment_info.l_seq + 1) // 2))
            if alignment_info.l_seq > 0
            else b"*"
        )
        if alignment_info.l_seq % 2:
            # When length is odd, the last nibble is junk
            seq = seq[:-1]

        qual = block.read(alignment_info.l_seq) if alignment_info.l_seq > 0 else b"*"
        if qual.startswith(b"\xFF"):
            # If qual is missing, it's filled with just \xFF repeated alignment_info.l_seq times
            qual = b"*"

        alignment_ref_name = (
            self._references[alignment_info.ref_id].name
            if alignment_info.ref_id != -1
            else "*"
        )
        next_ref_name = (
            self._references[alignment_info.next_ref_id].name
            if alignment_info.next_ref_id != -1
            else "*"
        )
        # aux = block.read()

        sam_line = "\t".join(
            [
                read_name.decode(),
                str(alignment_info.flag),
                str(alignment_ref_name),
                str(alignment_info.pos + 1),
                str(alignment_info.mapq),
                cigar.decode(),
                next_ref_name,
                str(alignment_info.next_pos + 1),
                str(alignment_info.tlen),
                seq.translate(self.seq_encoding).decode(),
                qual.translate(self.qual_mapping).decode(),
            ]
        )
        return sam_line


class AlignmentIterator(sam.AlignmentIterator):
    """
    BAM file parser (currently utilizing the SAM parser)
    """

    fmt = "BAM"

    def __init__(self, source):
        """Source should be a binary source accepted by BgzfReader"""
        self.source = source
        self._stream = _BAMToSamReader(source)
        self._index = 0
        self._read_header(self._stream)
