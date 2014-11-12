# Copyright 2002 by Yves Bastide and Brad Chapman.
# Copyright 2007 by Sebastian Bassi
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Functions to calculate assorted sequence checksums."""

# crc32, crc64, gcg, and seguid
# crc64 is adapted from BioPerl

from __future__ import print_function

from binascii import crc32 as _crc32
from Bio._py3k import _as_bytes

__docformat__ = "restructuredtext en"


def crc32(seq):
    """Returns the crc32 checksum for a sequence (string or Seq object)."""
    # NOTE - On Python 2 returns a signed int, on Python 3 it is unsigned
    # Docs suggest should use crc32(x) & 0xffffffff for consistency.
    # TODO - Should we return crc32(x) & 0xffffffff here?
    try:
        # Assume its a Seq object
        return _crc32(_as_bytes(str(seq)))
    except AttributeError:
        # Assume its a string/unicode
        return _crc32(_as_bytes(seq))


def _init_table_h():
    _table_h = []
    for i in range(256):
        l = i
        part_h = 0
        for j in range(8):
            rflag = l & 1
            l >>= 1
            if part_h & 1:
                l |= (1 << 31)
            part_h >>= 1
            if rflag:
                part_h ^= 0xd8000000
        _table_h.append(part_h)
    return _table_h

# Initialisation
_table_h = _init_table_h()


def crc64(s):
    """Returns the crc64 checksum for a sequence (string or Seq object)."""
    crcl = 0
    crch = 0
    for c in s:
        shr = (crch & 0xFF) << 24
        temp1h = crch >> 8
        temp1l = (crcl >> 8) | shr
        idx = (crcl ^ ord(c)) & 0xFF
        crch = temp1h ^ _table_h[idx]
        crcl = temp1l

    return "CRC-%08X%08X" % (crch, crcl)


def gcg(seq):
    """Returns the GCG checksum (int) for a sequence (string or Seq object).

    Given a nucleotide or amino-acid secuence (or any string),
    returns the GCG checksum (int). Checksum used by GCG program.
    seq type = str.

    Based on BioPerl GCG_checksum. Adapted by Sebastian Bassi
    with the help of John Lenton, Pablo Ziliani, and Gabriel Genellina.

    All sequences are converted to uppercase.
    """
    try:
        # Assume its a Seq object
        seq = str(seq)
    except AttributeError:
        # Assume its a string
        pass
    index = checksum = 0
    for char in seq:
        index += 1
        checksum += index * ord(char.upper())
        if index == 57:
            index = 0
    return checksum % 10000


def seguid(seq):
    """Returns the SEGUID (string) for a sequence (string or Seq object).

    Given a nucleotide or amino-acid secuence (or any string),
    returns the SEGUID string (A SEquence Globally Unique IDentifier).
    seq type = str.

    For more information about SEGUID, see:
    http://bioinformatics.anl.gov/seguid/
    DOI: 10.1002/pmic.200600032
    """
    import hashlib
    import base64
    m = hashlib.sha1()
    try:
        # Assume it's a Seq object
        seq = str(seq)
    except AttributeError:
        # Assume it's a string
        pass
    m.update(_as_bytes(seq.upper()))
    try:
        # For Python 3+
        return base64.encodebytes(m.digest()).decode().replace("\n", "").rstrip("=")
    except AttributeError:
        pass
    # For all other Pythons
    return base64.b64encode(m.digest()).rstrip("=")


if __name__ == "__main__":
    print("Quick self test")

    str_light_chain_one = "QSALTQPASVSGSPGQSITISCTGTSSDVGSYNLVSWYQQHPGK" \
                    + "APKLMIYEGSKRPSGVSNRFSGSKSGNTASLTISGLQAEDEADY" \
                    + "YCSSYAGSSTLVFGGGTKLTVL"

    str_light_chain_two = "QSALTQPASVSGSPGQSITISCTGTSSDVGSYNLVSWYQQHPGK" \
                    + "APKLMIYEGSKRPSGVSNRFSGSKSGNTASLTISGLQAEDEADY" \
                    + "YCCSYAGSSTWVFGGGTKLTVL"

    assert crc64(str_light_chain_one) == crc64(str_light_chain_two)
    assert 'CRC-44CAAD88706CC153' == crc64(str_light_chain_one)

    assert 'BpBeDdcNUYNsdk46JoJdw7Pd3BI' == seguid(str_light_chain_one)
    assert 'X5XEaayob1nZLOc7eVT9qyczarY' == seguid(str_light_chain_two)

    print("Done")
