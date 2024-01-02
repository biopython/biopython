"""
This module defines the dna-2 encoding scheme.

The dna-2 encoding scheme is experimental and subject to change or removal.

Do not import this module directly. Instead, import Bio.codecs.

Classes:
 - Codec - Defines the dna-2 encoding scheme.

Usage:
>>> import Bio.codecs
>>> "GATGCGATCGTAGCTGAT".encode("dna-2").decode("dna-2")
'GATGCGATCGTAGCTGAT'
"""

import Bio._utils
import codecs
import math

from typing import Union

_nucleotide_characters = set("ATCGatcg")


def _count_nucleotide_characters(string: str, should_error: bool) -> int:
    nucleotide_count = 0

    for character in string:
        if character in _nucleotide_characters:
            nucleotide_count += 1
        elif should_error:
            raise ValueError(f"Unexpected character '{character}'")

    return nucleotide_count


def _encode_nucleotide(nucleotide: str) -> int:
    assert len(nucleotide) == 1
    assert nucleotide in _nucleotide_characters

    if nucleotide in "Tt":
        return 0b00
    elif nucleotide in "Cc":
        return 0b01
    elif nucleotide in "Aa":
        return 0b10
    else:
        assert nucleotide in "Gg"
        return 0b11


def _decode_nucleotide(byte: int) -> str:
    return {
        0b00: "T",
        0b01: "C",
        0b10: "A",
        0b11: "G",
    }[byte]


def _encode_nucleotides_into_byte(nucleotides: Union[list[str], str]) -> int:
    nucleotide_count = len(nucleotides)
    assert nucleotide_count <= 4
    for nucleotide in nucleotides:
        assert nucleotide in _nucleotide_characters

    bit_shifts = [6, 4, 2, 0][0:nucleotide_count]
    result = 0

    for nucleotide_index, bit_shift in enumerate(bit_shifts):
        nucleotide = nucleotides[nucleotide_index]
        result |= _encode_nucleotide(nucleotide) << bit_shift

    if nucleotide_count < 4:
        result |= nucleotide_count

    assert 0x0 <= result <= 0xFF

    return result


def _decode_nucleotides_from_byte(byte: int, is_last_byte: bool = False) -> str:
    assert 0x0 <= byte <= 0xFF
    nucleotide_count = byte & 0b11 if is_last_byte else 4
    assert 0 <= nucleotide_count <= 4

    nucleotides = []
    bit_shifts = [6, 4, 2, 0][0:nucleotide_count]

    for bit_shift in bit_shifts:
        nucleotide = _decode_nucleotide(byte >> bit_shift & 0b11)
        nucleotides.append(nucleotide)

    return "".join(nucleotides)


class Codec(codecs.Codec):
    """Defines the dna-2 encoding scheme."""

    def encode(self, input, errors="strict"):
        """Encode the input string representing a DNA sequence into binary using the dna-2 encoding scheme."""
        assert isinstance(input, str)
        assert errors in {"strict", "ignore"}

        nucleotide_count = _count_nucleotide_characters(
            input, should_error=errors == "strict"
        )
        bit_count = 2 * (nucleotide_count + 1)
        byte_count: int = math.ceil(bit_count / 8)

        result = bytearray(byte_count)
        byte_index = 0
        nucleotides = []

        for character in input:
            assert len(nucleotides) <= 4
            if character not in _nucleotide_characters:
                continue

            nucleotide = character
            nucleotides.append(nucleotide)

            if len(nucleotides) == 4:
                result[byte_index] = _encode_nucleotides_into_byte(nucleotides)
                byte_index += 1
                nucleotides = []

        assert byte_index == byte_count - 1
        result[byte_index] = _encode_nucleotides_into_byte(nucleotides)

        return bytes(result), len(input)

    def decode(self, input, errors="strict"):
        """Decode the input bytes into a string according to the dna-2 encoding scheme."""
        assert (
            isinstance(input, bytes)
            or isinstance(input, bytearray)
            or isinstance(input, memoryview)
        )

        strings = [_decode_nucleotides_from_byte(byte) for byte in input[0:-1]] + [
            _decode_nucleotides_from_byte(input[-1], is_last_byte=True)
        ]

        return "".join(strings), len(input)


if __name__ == "__main__":
    Bio._utils.run_doctest()
