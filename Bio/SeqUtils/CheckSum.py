# crc32, crc64, gcg, and seguid
# crc64 is adapted from BioPerl

from binascii import crc32

def _init_table_h():
    _table_h = []
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
    return _table_h

# Initialisation
_table_h = _init_table_h()

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


def gcg(seq):
    """ given a nucleotide or amino-acid secuence (or any string),
        returns the GCG checksum (int). Checksum used by GCG program.
        seq type = str.
        Based on BioPerl GCG_checksum. Adapted by Sebastian Bassi
        with the help of John Lenton and Pablo Ziliani.
        Python 2.5 Oneliner by Gabriel Genellina.
        All sequences are converted to uppercase """
    import sys
    if sys.version_info >= (2,4):
        # GCG_checksum version for Python 2.4 and up.
        # By Sebastian Bassi and Gabriel Genellina.
        from itertools import cycle, izip
        return sum(n*ord(c.upper()) for (n,c) in izip(cycle(range(1,58)),seq)) % 10000
    else:
        #slower version for Python 2.3
        index = checksum = 0
        if type(seq)!=type("aa"):
            seq=seq.tostring().upper()
        else:
            seq=seq.upper()
        for char in seq.upper():
            index += 1
            checksum += index * ord(char)
            if index == 57: index = 0
        return checksum % 10000


def seguid(seq):
    """ given a nucleotide or amino-acid secuence (or any string),
        returns the SEGUID string (A SEquence Globally Unique IDentifier).
        seq type = str. 
        For more information about SEGUID, see:
        http://bioinformatics.anl.gov/seguid/
        DOI: 10.1002/pmic.200600032 """
    try:
        #Python 2.5 sha1 is in hashlib
        import hashlib
        m = hashlib.sha1()
    except:
        #For older versions 
        import sha
        m = sha.new()
    import base64
    if type(seq)!=type("aa"):
        seq=seq.tostring().upper()
    else:
        seq=seq.upper()
    m.update(seq)
    try:
        #For Python 2.5
        return base64.b64encode(m.digest()).rstrip("=")
    except:
        #For older versions
        import os
        return base64.encodestring(m.digest()).replace(os.linesep,"").rstrip("=")
