# Copyright 2003 by Bartek Wilczynski.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Parsing TRANSFAC files
"""

from Bio.Motif import Motif as BaseMotif

class Motif(BaseMotif, dict):
    """A Bio.Motif.TRANSFAC.Motif stores the information in one TRANSFAC
motif. This class inherits from the Bio.Motif.Motif base class, as well
as from a Python dictionary. All motif information found by the parser
is stored as attributes of the base class when possible; see the
Bio.Motif.Motif base class for a description of these attributes. All
other information associated with the motif is stored as (key, value)
pairs in the dictionary, where the key is the two-letter fields as found
in the TRANSFAC file.

These fields are commonly found in TRANSFAC files:
    AC:    Accession number
    AS:    Accession numbers, secondary
    BA:    Statistical basis
    BF:    Binding factors
    BS:    Factor binding sites underlying the matrix
           [SITE accession number; start position for matrix sequence;
            length of sequence used; number of gaps inserted; strand
            orientation]
    CC:    Comments
    CO:    Copyright notice
    DE:    Short factor description
    DR:    External databases
           [EMBL/GEnBank accession number; EMBL identifier (1st: last
            position of the TRANSFAC sequence element)]
    DT:    Date created/updated
    HC:    Subfamilies
    HP:    Superfamilies
    ID:    Identifier
    NA:    Name of the binding factor
    OC:    Taxonomic classification
    OS:    Species/Taxon
    OV:    Older version
    PV:    Preferred version
    RA:    Reference authors
    RL:    Reference data
    RN:    Reference number
    RT:    Reference title
    RX:    PubMed ID
    TY:    Type
    XX:    Empty line; these are not stored in the Record.

For more information, see the TRANSFAC documentation.
"""

class Record(object):
    """A Bio.Motif.TRANSFAC.Record stores the information in a TRANSFAC
matrix table.

Attributes:
    o version:   The version number, corresponding to the 'VV' field
                 in the TRANSFAC file;
    o motifs:    The list of motifs.
"""

def read(handle):
    """record = read(handle)"""
    record = Record()
    for line in handle:
        line = line.strip()
        if line=='//':
            break
        elif line=='XX':
            pass
        else:
            key, value = line.split(None, 1)
            assert key=='VV'
            record.version = value
    record.motifs = []
    motif = None
    status = None
    for line in handle:
        line = line.strip()
        if motif==None:
            motif = Motif()
        if line=='//':
            record.motifs.append(motif)
            motif = None
            status = None
        elif line=='XX':
            pass
        else:
            key, value = line.split(None, 1)
            if status=="freq":
               try:
                   i = int(key)
               except ValueError:
                   status = None
               else:
                   motif.length+=1
                   assert i==motif.length
                   values = value.split()
                   for c, v in zip("ACGT", values):
                       motif.counts[c].append(float(v))
                   continue
            if key=='P0':
                assert status!="freq"
                assert motif.counts==None
                motif.counts = {}
                assert value.split()[:4]==['A','C','G','T']
                motif.length = 0
                for c in "ACGT":
                    motif.counts[c] = []
                status = "freq"
            else:
                motif[key] = value
    return record
