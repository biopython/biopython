# Copyright 2008 by Bartek Wilczynski
# Adapted from  Bio.MEME.Parser by Jason A. Hackney.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

from __future__ import print_function

from Bio.Alphabet import IUPAC
from Bio import Seq
from Bio import motifs


def read(handle):
    """Parses the text output of the MEME program into a meme.Record object.

    Example:

    >>> from Bio.motifs import meme
    >>> with open("meme.output.txt") as f:
    ...     record = meme.read(f)
    >>> for motif in record:
    ...     for instance in motif.instances:
    ...         print(instance.motif_name, instance.sequence_name, instance.strand, instance.pvalue)

    """
    record = Record()
    __read_version(record, handle)
    __read_datafile(record, handle)
    __read_alphabet(record, handle)
    __read_sequences(record, handle)
    __read_command(record, handle)
    for line in handle:
        if line.startswith('MOTIF  1'):
            break
    else:
        raise ValueError('Unexpected end of stream')
    alphabet = record.alphabet
    revcomp = 'revcomp' in record.command
    while True:
        length, num_occurrences, evalue = __read_motif_statistics(line)
        name = __read_motif_name(handle)
        instances = __read_motif_sequences(handle, name, alphabet, length, revcomp)
        motif = Motif(alphabet, instances)
        motif.length = length
        motif.num_occurrences = num_occurrences
        motif.evalue = evalue
        motif.name = name
        record.append(motif)
        __skip_unused_lines(handle)
        try:
            line = next(handle)
        except StopIteration:
            raise ValueError('Unexpected end of stream: Expected to find new motif, or the summary of motifs')
        if line.startswith("SUMMARY OF MOTIFS"):
            break
        if not line.startswith('MOTIF'):
            raise ValueError("Line does not start with 'MOTIF':\n%s" % line)
    return record


class Motif(motifs.Motif):
    """A subclass of Motif used in parsing MEME (and MAST) output.

    This subclass defines functions and data specific to MEME motifs.
    This includes the motif name, the evalue for a motif, and its number
    of occurrences.
    """
    def __init__(self, alphabet=None, instances=None):
        motifs.Motif.__init__(self, alphabet, instances)
        self.evalue = 0.0
        self.num_occurrences = 0
        self.name = None


class Instance(Seq.Seq):
    """A class describing the instances of a MEME motif, and the data thereof.
    """
    def __init__(self, *args, **kwds):
        Seq.Seq.__init__(self, *args, **kwds)
        self.sequence_name = ""
        self.start = 0
        self.pvalue = 1.0
        self.strand = 0
        self.length = 0
        self.motif_name = ""


class Record(list):
    """A class for holding the results of a MEME run.

    A meme.Record is an object that holds the results from running
    MEME. It implements no methods of its own.

    The meme.Record class inherits from list, so you can access individual
    motifs in the record by their index. Alternatively, you can find a motif
    by its name:

    >>> from Bio import motifs
    >>> with open("meme.output.txt") as f:
    ...     record = motifs.parse(f, 'MEME')
    >>> motif = record[0]
    >>> print(motif.name)
    Motif 1
    >>> motif = record['Motif 1']
    >>> print(motif.name)
    Motif 1
    """

    def __init__(self):
        """__init__ (self)"""
        self.version = ""
        self.datafile = ""
        self.command = ""
        self.alphabet = None
        self.sequences = []

    def __getitem__(self, key):
        if isinstance(key, str):
            for motif in self:
                if motif.name==key:
                    return motif
        else:
            return list.__getitem__(self, key)


# Everything below is private


def __read_version(record, handle):
    for line in handle:
        if line.startswith('MEME version'):
            break
    else:
        raise ValueError("Improper input file. File should contain a line starting MEME version.")
    line = line.strip()
    ls = line.split()
    record.version = ls[2]


def __read_datafile(record, handle):
    for line in handle:
        if line.startswith('TRAINING SET'):
            break
    else:
        raise ValueError("Unexpected end of stream: 'TRAINING SET' not found.")
    try:
        line = next(handle)
    except StopIteration:
        raise ValueError("Unexpected end of stream: Expected to find line starting with '****'")
    if not line.startswith('****'):
        raise ValueError("Line does not start with '****':\n%s" % line)
    try:
        line = next(handle)
    except StopIteration:
        raise ValueError("Unexpected end of stream: Expected to find line starting with 'DATAFILE'")
    if not line.startswith('DATAFILE'):
        raise ValueError("Line does not start with 'DATAFILE':\n%s" % line)
    line = line.strip()
    line = line.replace('DATAFILE= ', '')
    record.datafile = line


def __read_alphabet(record, handle):
    try:
        line = next(handle)
    except StopIteration:
        raise ValueError("Unexpected end of stream: Expected to find line starting with 'ALPHABET'")
    if not line.startswith('ALPHABET'):
        raise ValueError("Line does not start with 'ALPHABET':\n%s" % line)
    line = line.strip()
    line = line.replace('ALPHABET= ', '')
    if line == 'ACGT':
        al = IUPAC.unambiguous_dna
    else:
        al = IUPAC.protein
    record.alphabet = al


def __read_sequences(record, handle):
    try:
        line = next(handle)
    except StopIteration:
        raise ValueError("Unexpected end of stream: Expected to find line starting with 'Sequence name'")
    if not line.startswith('Sequence name'):
        raise ValueError("Line does not start with 'Sequence name':\n%s" % line)
    try:
        line = next(handle)
    except StopIteration:
        raise ValueError("Unexpected end of stream: Expected to find line starting with '----'")
    if not line.startswith('----'):
        raise ValueError("Line does not start with '----':\n%s" % line)
    for line in handle:
        if line.startswith('***'):
            break
        line = line.strip()
        ls = line.split()
        record.sequences.append(ls[0])
        if len(ls) == 6:
            record.sequences.append(ls[3])
    else:
        raise ValueError("Unexpected end of stream: Expected to find line starting with '***'")


def __read_command(record, handle):
    for line in handle:
        if line.startswith('command:'):
            break
    else:
        raise ValueError("Unexpected end of stream: Expected to find line starting with 'command'")
    line = line.strip()
    line = line.replace('command: ', '')
    record.command = line


def __read_motif_statistics(line):
    line = line[5:].strip()
    ls = line.split()
    length = int(ls[3])
    num_occurrences = int(ls[6])
    evalue = float(ls[12])
    return length, num_occurrences, evalue


def __read_motif_name(handle):
    for line in handle:
        if 'sorted by position p-value' in line:
            break
    else:
        raise ValueError('Unexpected end of stream: Failed to find motif name')
    line = line.strip()
    words = line.split()
    name = " ".join(words[0:2])
    return name


def __read_motif_sequences(handle, motif_name, alphabet, length, revcomp):
    try:
        line = next(handle)
    except StopIteration:
        raise ValueError('Unexpected end of stream: Failed to find motif sequences')
    if not line.startswith('---'):
        raise ValueError("Line does not start with '---':\n%s" % line)
    try:
        line = next(handle)
    except StopIteration:
        raise ValueError("Unexpected end of stream: Expected to find line starting with 'Sequence name'")
    if not line.startswith('Sequence name'):
        raise ValueError("Line does not start with 'Sequence name':\n%s" % line)
    try:
        line = next(handle)
    except StopIteration:
        raise ValueError('Unexpected end of stream: Failed to find motif sequences')
    if not line.startswith('---'):
        raise ValueError("Line does not start with '---':\n%s" % line)
    instances = []
    for line in handle:
        if line.startswith('---'):
            break
        line = line.strip()
        words = line.split()
        if revcomp:
            strand = words.pop(1)
        else:
            strand = '+'
        sequence = words[4]
        assert len(sequence)==length
        instance = Instance(sequence, alphabet)
        instance.motif_name = motif_name
        instance.sequence_name = words[0]
        instance.start = int(words[1])
        instance.pvalue = float(words[2])
        instance.strand = strand
        instance.length = length
        instances.append(instance)
    else:
        raise ValueError('Unexpected end of stream')
    return motifs.Instances(instances, alphabet)


def __skip_unused_lines(handle):
    for line in handle:
        if line.startswith('log-odds matrix'):
            break
    else:
        raise ValueError("Unexpected end of stream: Expected to find line starting with 'log-odds matrix'")
    for line in handle:
        if line.startswith('---'):
            break
    else:
        raise ValueError("Unexpected end of stream: Expected to find line starting with '---'")
    for line in handle:
        if line.startswith('letter-probability matrix'):
            break
    else:
        raise ValueError("Unexpected end of stream: Expected to find line starting with 'letter-probability matrix'")
    for line in handle:
        if line.startswith('---'):
            break
    else:
        raise ValueError("Unexpected end of stream: Expected to find line starting with '---'")
    for line in handle:
        if line.startswith('Time'):
            break
    else:
        raise ValueError("Unexpected end of stream: Expected to find line starting with 'Time'")
    try:
        line = next(handle)
    except StopIteration:
        raise ValueError('Unexpected end of stream: Expected to find blank line')
    if line.strip():
        raise ValueError("Expected blank line, but got:\n%s" % line)
    try:
        line = next(handle)
    except StopIteration:
        raise ValueError("Unexpected end of stream: Expected to find line starting with '***'")
    if not line.startswith('***'):
        raise ValueError("Line does not start with '***':\n%s" % line)
    for line in handle:
        if line.strip():
            break
    else:
        raise ValueError("Unexpected end of stream: Expected to find line starting with '***'")
    if not line.startswith('***'):
        raise ValueError("Line does not start with '***':\n%s" % line)
