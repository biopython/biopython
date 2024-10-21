"""Bio.Align support for the "bam" pairwise alignment format."""

import numpy as np
from Bio import bgzf
from io import StringIO, BytesIO
from typing import NamedTuple, Optional, Literal, Iterator, Any, Union
import struct
import array
from Bio.Align import sam
import os


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


class _TagData(NamedTuple):
    type: bytes
    data: Union[bytes, int, float, array.array]


class _BAMToSamReader(bgzf.BgzfReader):
    """Limited translation reader from BAM to SAM which allows seek(0), read(0), and iteration"""

    _references: Optional[list[_ReferenceMetadata]]
    _header_iterator: Optional[Iterator[str]]

    seq_encoding: bytes = bytes.maketrans(bytes(range(16)), b"=ACMGRSVTWYHKDBN")
    cigar_encoding: bytes = bytes.maketrans(b"012345678", b"MIDNSHP=X")
    read_lengths: dict[bytes, int] = {
        b"A": 1,
        b"c": 1,
        b"C": 1,
        b"s": 2,
        b"S": 2,
        b"i": 4,
        b"I": 4,
        b"f": 4,
    }

    type_format_mapping: dict[bytes, bytes] = {
        b"A": b"c",
        b"c": b"<b",
        b"C": b"<B",
        b"s": b"<h",
        b"S": b"<H",
        b"i": b"<i",
        b"I": b"<I",
        b"f": b"f",
    }
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

    @classmethod
    def parse_aux_data(cls, aux_bytes: BytesIO) -> dict[bytes, _TagData]:

        aux_data = {}
        tag = aux_bytes.read(2)
        while tag:
            # Ends when aux_bytes is exhausted, so aux_bytes.read(2) == b""
            data_type = aux_bytes.read(1)

            read_length = cls.read_lengths.get(data_type)
            if read_length is not None:
                # Single numeric type, like tag|type|data
                python_data_type = cls.type_format_mapping[data_type]
                data_bytes = aux_bytes.read(read_length)
                data = struct.unpack(python_data_type, data_bytes)[0]
                aux_data[tag] = _TagData(data_type, data)
            elif data_type in [b"Z", b"H"]:
                # Null terminated character or hex data. We're relying on the fact that aux_bytes is a BytesIO object,
                # and so the buffer is everything. Otherwise, we would need to slowly look byte by byte
                data_start_index = aux_bytes.tell()
                underlying_data = aux_bytes.getbuffer().tobytes()
                data_end_index = underlying_data.find(b"\00", data_start_index)
                tag_data = (
                    aux_bytes.read()
                    if data_end_index == -1
                    else aux_bytes.read(data_end_index - data_start_index)
                )

                aux_data[tag] = _TagData(data_type, tag_data)
                # Jump forward by one to ignore terminating null byte
                aux_bytes.seek(1, os.SEEK_CUR)
            elif data_type == b"B":
                # Array type, composed like: tag|b'B'|entry_type|length|data
                entry_type = aux_bytes.read(1)
                python_entry_type = cls.type_format_mapping[entry_type].decode()

                entry_size = cls.read_lengths[entry_type]
                array_len_bytes = aux_bytes.read(4)
                array_len = _to_uint(array_len_bytes)

                array_data = array.array(
                    python_entry_type, aux_bytes.read(array_len * entry_size)
                )
                aux_data[tag] = _TagData(data_type + entry_type, array_data)
            else:
                raise ValueError(
                    f"Unexpected tag type {data_type.decode()} for tag [{tag.decode()}]"
                )

            tag = aux_bytes.read(2)

        return aux_data

    @classmethod
    def translate_aux_data(cls, tag: bytes, data_type: bytes, data: bytes) -> str:
        if data_type in (b"A", b"Z", b"H"):
            data_string = data.decode()
        elif data_type[0:1] == b"B":
            entry_type = data_type[1:2]
            entry_type = b"f" if entry_type == b"f" else b"i"
            data_string = (entry_type + b",".join(np.char.array(data))).decode()
        elif data_type in (b"c", b"C", b"s", b"S", b"i", b"I", b"f"):
            data_type = b"f" if data_type == b"f" else b"i"
            data_string = str(data)
        else:
            raise ValueError
        return f"{tag.decode()}:{data_type.decode()}:{data_string}"

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

        reference_metadata = self._references[alignment_info.ref_id]
        alignment_ref_name = (
            reference_metadata.name if alignment_info.ref_id != -1 else "*"
        )
        next_ref_name = (
            self._references[alignment_info.next_ref_id].name
            if alignment_info.next_ref_id != -1
            else "*"
        )
        aux_data = self.parse_aux_data(block)

        placeholder = (
            f"{alignment_info.l_seq}S{reference_metadata.length}N".encode()
            if alignment_info.ref_id != -1
            else None
        )
        if cigar == placeholder and b"CG" in aux_data:
            # Check if cigar operations are in the tag data which happens if there are more than 65535 CIGAR operations
            # cigar = self.get_tag_data(aux, b"CG", b"BI")
            cigar = self.cigar_nibble_decode(aux_data[b"CG"].data.tobytes())
            del aux_data[b"CG"]

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
            + [
                self.translate_aux_data(tag, data_type, data)
                for tag, (data_type, data) in aux_data.items()
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
