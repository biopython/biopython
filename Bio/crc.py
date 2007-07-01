# crc32 and crc64
# crc64 is adapted from BioPerl

import warnings
warnings.warn("Bio.crc is deprecated; use crc32 and crc64 in Bio.SeqUtils.CheckSum instead", DeprecationWarning)


from binascii import crc32

_table_h = []

def crc64(s):
    crcl = 0
    crch = 0
    for c in s:
        shr = (crch & 0xFF) << 24
        temp1h = crch >> 8
        temp1l = (crcl >> 8) | shr
        idx  = (crcl ^ ord(c)) & 0xFF
        crch = temp1h ^ _table_h[idx]
        crcl = temp1l

    return "CRC-%08X%08X" % (crch, crcl)

# Initialisation
for i in range(256):
    l = i
    part_h = 0
    for j in range(8):
        rflag = l & 1
        l >>= 1
        if part_h & 1: l |= (1L << 31)
        part_h >>= 1L
        if rflag: part_h ^= 0xd8000000L
    _table_h.append(part_h)

