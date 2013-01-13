# Copyright 2003 by Bartek Wilczynski.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Parsing TRANSFAC files
"""

import warnings
from Bio import BiopythonExperimentalWarning

warnings.warn("Bio.Motif.TRANSFAC is experimental code. While it is usable, \
    the code is subject to change without warning", BiopythonExperimentalWarning)

from Bio.Motif import NewMotif as BaseMotif
from Bio.Alphabet import IUPAC


class Motif(BaseMotif, dict):
    """A Bio.Motif.TRANSFAC.Motif stores the information in one TRANSFAC
motif. This class inherits from the Bio.Motif.Motif base class, as well
as from a Python dictionary. All motif information found by the parser
is stored as attributes of the base class when possible; see the
Bio.Motif.Motif base class for a description of these attributes. All
other information associated with the motif is stored as (key, value)
pairs in the dictionary, where the key is the two-letter fields as found
in the TRANSFAC file. References are an exception: These are stored in
the .references attribute.

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
    TY:    Type
    XX:    Empty line; these are not stored in the Record.

References are stored in an .references attribute, which is a list of
dictionaries with the following keys:
    RN:    Reference number
    RA:    Reference authors
    RL:    Reference data
    RT:    Reference title
    RX:    PubMed ID

For more information, see the TRANSFAC documentation.
"""
    multiple_value_keys = set(['BF', 'OV', 'HP', 'BS', 'HC', 'DT', 'DR'])
    # These keys can occur multiple times for one motif

    reference_keys = set(['RX', 'RA', 'RT', 'RL'])
    # These keys occur for references


class Record(list):
    """A Bio.Motif.TRANSFAC.Record stores the information in a TRANSFAC
matrix table. The record inherits from a list containing the individual
motifs.

Attributes:
    o version:   The version number, corresponding to the 'VV' field
                 in the TRANSFAC file;
"""
    def __init__(self):
        self.version = None

    @property
    def motifs(self):
        import warnings
        warnings.warn("""\
The .motifs attribute is now obsolete, and will be deprecated and removed
in a future release of Biopython. This class now inherits from list, so
instead of record.motifs[i], please use record[i].
""", PendingDeprecationWarning)
        return self

    def __str__(self):
        blocks = []
        if self.version is not None:
            block = """\
VV  %s
XX
//
""" % self.version
            blocks.append(block)
        for motif in self:
            block = str(motif)
            blocks.append(block)
        text = "".join(blocks)
        return text


def read(handle):
    """record = read(handle)"""
    annotations = {}
    references = []
    counts = None
    record = Record()
    for line in handle:
        line = line.strip()
        key, value = line[:2], line[4:]
        if key=='VV':
            record.version = value
        elif key=='P0':
            counts = {}
            assert value.split()[:4]==['A','C','G','T']
            length = 0
            for c in "ACGT":
                counts[c] = []
            for line in handle:
                key, value = line[:2], line[4:]
                try:
                    i = int(key)
                except ValueError:
                    break
                length+=1
                assert i==length
                values = value.split()
                for c, v in zip("ACGT", values):
                    counts[c].append(float(v))
        if line=='XX':
            pass
        elif key=='RN':
            index, accession = value.split(";")
            assert index[0]=='['
            assert index[-1]==']'
            index = int(index[1:-1])
            assert len(references)==index-1
            reference = {key: value}
            references.append(reference)
        elif key=='//':
            if counts is not None:
                motif = Motif(alphabet=IUPAC.unambiguous_dna, counts=counts)
                motif.update(annotations)
                motif.references = references
                record.append(motif)
            annotations = {}
            references = []
        elif key in Motif.reference_keys:
            reference[key] = value
        elif key in Motif.multiple_value_keys:
            if not key in annotations:
                annotations[key] = []
            annotations[key].append(value)
        else:
            annotations[key] = value
    return record
