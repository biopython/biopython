# Copyright 2002 by Yves Bastide and Brad Chapman.
# Copyright 2007 by Sebastian Bassi
# All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Functions to calculate assorted sequence checksums."""

# crc32, crc64, gcg, and seguid
# crc64 is adapted from BioPerl


import binascii


def crc32(seq):
    """Return the crc32 checksum for a sequence (string or Seq object).

    Note that the case is important:

    >>> crc32("ACGTACGTACGT")
    20049947
    >>> crc32("acgtACGTacgt")
    1688586483

    """
    try:
        # Assume it's a Seq object
        s = bytes(seq)
    except TypeError:
        # Assume it's a string
        s = seq.encode()
    return binascii.crc32(s)


def _init_table_h():
    _table_h = []
    for i in range(256):
        part_l = i
        part_h = 0
        for j in range(8):
            rflag = part_l & 1
            part_l >>= 1
            if part_h & 1:
                part_l |= 1 << 31
            part_h >>= 1
            if rflag:
                part_h ^= 0xD8000000
        _table_h.append(part_h)
    return _table_h


# Initialisation
_table_h = _init_table_h()


def crc64(s):
    """Return the crc64 checksum for a sequence (string or Seq object).

    Note that the case is important:

    >>> crc64("ACGTACGTACGT")
    'CRC-C4FBB762C4A87EBD'
    >>> crc64("acgtACGTacgt")
    'CRC-DA4509DC64A87EBD'

    """
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
    """Return the GCG checksum (int) for a sequence (string or Seq object).

    Given a nucleotide or amino-acid sequence (or any string),
    returns the GCG checksum (int). Checksum used by GCG program.
    seq type = str.

    Based on BioPerl GCG_checksum. Adapted by Sebastian Bassi
    with the help of John Lenton, Pablo Ziliani, and Gabriel Genellina.

    All sequences are converted to uppercase.

    >>> gcg("ACGTACGTACGT")
    5688
    >>> gcg("acgtACGTacgt")
    5688

    """
    index = checksum = 0
    for char in seq:
        index += 1
        checksum += index * ord(char.upper())
        if index == 57:
            index = 0
    return checksum % 10000


def seguid(seq):
    """Return the SEGUID (string) for a sequence (string or Seq object).

    Given a nucleotide or amino-acid sequence (or any string),
    returns the SEGUID string (A SEquence Globally Unique IDentifier).
    seq type = str.

    Note that the case is not important:

    >>> seguid("ACGTACGTACGT")
    'If6HIvcnRSQDVNiAoefAzySc6i4'
    >>> seguid("acgtACGTacgt")
    'If6HIvcnRSQDVNiAoefAzySc6i4'

    For more information about SEGUID, see:
    http://bioinformatics.anl.gov/seguid/
    https://doi.org/10.1002/pmic.200600032
    """
    import hashlib
    import base64

    m = hashlib.sha1()
    try:
        # Assume it's a Seq object
        seq = bytes(seq)
    except TypeError:
        # Assume it's a string
        seq = seq.encode()
    m.update(seq.upper())
    tmp = base64.encodebytes(m.digest())
    return tmp.decode().replace("\n", "").rstrip("=")


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
