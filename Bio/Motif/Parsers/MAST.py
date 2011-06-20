# Copyright 2008 by Bartek Wilczynski.
# Adapted from Bio.MEME.Parser by Jason A. Hackney.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from Bio.Alphabet import IUPAC
from Bio.Motif.Parsers.MEME import MEMEMotif


class Record(object):
    """The class for holding the results from a MAST run.
    
    A MAST.Record holds data about matches between motifs and sequences.
    The motifs held by the Record are objects of the class MEMEMotif.
    
    Methods:
    get_motif_by_name (motif_name): returns a MEMEMotif with the given
    name.
    """

    def __init__ (self):
        self.sequences = []
        self.version = ""
        self.database = ""
        self.diagrams = {}
        self.alphabet = None
        self.motifs = []
    
    def get_motif_by_name (self, name):
        for m in self.motifs:
            if m.name == name:
                return m

def read(handle):
    """read(handle)"""
    record = Record()
    __read_version(record, handle)
    __read_database_and_motifs(record, handle)
    __read_section_i(record, handle)
    __read_section_ii(record, handle)
    __read_section_iii(record, handle)
    return record


# Everything below is private


def __read_version(record, handle):
    for line in handle:
        if "MAST version" in line:
            break
    else:
        raise ValueError("Improper input file. Does not begin with a line with 'MAST version'")
    record.version = line.strip().split()[2]


def __read_database_and_motifs(record, handle):
    for line in handle:
        if line.startswith('DATABASE AND MOTIFS'):
            break
    line = handle.next()
    if not line.startswith('****'):
        raise ValueError("Line does not start with '****':\n%s" % line)
    line = handle.next()
    if not 'DATABASE' in line:
        raise ValueError("Line does not contain 'DATABASE':\n%s" % line)
    words = line.strip().split()
    record.database = words[1]
    if words[2] == '(nucleotide)':
        record.alphabet = IUPAC.unambiguous_dna
    elif words[2] == '(peptide)':
        record.alphabet = IUPAC.protein
    for line in handle:
        if 'MOTIF WIDTH' in line:
            break
    line = handle.next()
    if not '----' in line:
        raise ValueError("Line does not contain '----':\n%s" % line)
    for line in handle:
        if not line.strip():
            break
        words = line.strip().split()
        motif = MEMEMotif()
        motif.alphabet = record.alphabet
        motif.name = words[0]
        motif.length = int(words[1])
        # motif.add_instance(words[2])
        record.motifs.append(motif)


def __read_section_i(record, handle):
    for line in handle:
        if line.startswith('SECTION I:'):
            break
    for line in handle:
        if line.startswith('SEQUENCE NAME'):
            break
    line = handle.next()
    if not line.startswith('---'):
        raise ValueError("Line does not start with '---':\n%s" % line)
    for line in handle:
        if not line.strip():
            break
        else:
            sequence, description_evalue_length = line.split(None, 1)
            record.sequences.append(sequence)
    line = handle.next()
    if not line.startswith('****'):
        raise ValueError("Line does not start with '****':\n%s" % line)


def __read_section_ii(record, handle):
    for line in handle:
        if line.startswith('SECTION II:'):
            break
    for line in handle:
        if line.startswith('SEQUENCE NAME'):
            break
    line = handle.next()
    if not line.startswith('---'):
        raise ValueError("Line does not start with '---':\n%s" % line)
    for line in handle:
        if not line.strip():
            break
        elif line.startswith(" "):
            diagram = line.strip()
            record.diagrams[sequence] += diagram
        else:
            sequence, pvalue, diagram = line.split()
            record.diagrams[sequence] = diagram
    line = handle.next()
    if not line.startswith('****'):
        raise ValueError("Line does not start with '****':\n%s" % line)


def __read_section_iii(record, handle):
    for line in handle:
        if line.startswith('SECTION III:'):
            break
    for line in handle:
        if line.startswith('****'):
            break
    for line in handle:
        if line.startswith('*****'):
            break
    for line in handle:
        if line.strip():
            break
