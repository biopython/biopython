# Copyright 2008 by Bartek Wilczynski
# Adapted from  Bio.MEME.Parser by Jason A. Hackney.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Module for the support of MEME motif format."""

from __future__ import print_function

from Bio import Seq
from Bio import motifs


def read(handle):
    """Parse the text output of the MEME program into a meme.Record object.

    Examples
    --------
    >>> from Bio.motifs import meme
    >>> with open("motifs/meme.out") as f:
    ...     record = meme.read(f)
    >>> for motif in record:
    ...     for instance in motif.instances:
    ...         print(instance.motif_name, instance.sequence_name, instance.strand, instance.pvalue)
    Motif 1 SEQ10; + 8.71e-07
    Motif 1 SEQ9; + 8.71e-07
    Motif 1 SEQ8; + 8.71e-07
    Motif 1 SEQ7; + 8.71e-07
    Motif 1 SEQ6; + 8.71e-07
    Motif 1 SEQ5; + 8.71e-07
    Motif 1 SEQ4; + 8.71e-07
    Motif 1 SEQ3; + 8.71e-07
    Motif 1 SEQ2; + 8.71e-07
    Motif 1 SEQ1; + 8.71e-07

    """
    record = Record()
    __read_version(record, handle)
    __read_datafile(record, handle)
    __read_alphabet(record, handle)
    __read_sequences(record, handle)
    __read_command(record, handle)
    for line in handle:
        if line.startswith('MOTIF '):
            break
    else:
        raise ValueError('Unexpected end of stream')
    alphabet = record.alphabet
    revcomp = 'revcomp' in record.command
    while True:
        motif_number, length, num_occurrences, evalue = __read_motif_statistics(line)
        name = __read_motif_name(handle)
        instances = __read_motif_sequences(handle, name, alphabet, length, revcomp)
        motif = Motif(alphabet, instances)
        motif.length = length
        motif.num_occurrences = num_occurrences
        motif.evalue = evalue
        motif.name = name
        record.append(motif)
        assert len(record) == motif_number
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
        """Initialize the class."""
        motifs.Motif.__init__(self, alphabet, instances)
        self.evalue = 0.0
        self.num_occurrences = 0
        self.name = None
        self.id = None
        self.alt_id = None


class Instance(Seq.Seq):
    """A class describing the instances of a MEME motif, and the data thereof."""

    def __init__(self, *args, **kwds):
        """Initialize the class."""
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
    >>> with open("motifs/meme.out") as f:
    ...     record = motifs.parse(f, 'MEME')
    >>> motif = record[0]
    >>> print(motif.name)
    Motif 1
    >>> motif = record['Motif 1']
    >>> print(motif.name)
    Motif 1
    """

    def __init__(self):
        """Initialize."""
        self.version = ""
        self.datafile = ""
        self.command = ""
        self.alphabet = None
        self.sequences = []

    def __getitem__(self, key):
        """Return the motif of index key."""
        if isinstance(key, str):
            for motif in self:
                if motif.name == key:
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
        raise ValueError(
            "Unexpected end of stream: 'TRAINING SET' not found. This can happen with " +
            "minimal MEME files (MEME databases) which are not supported yet.")
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
        al = 'ACGT'
    elif line == 'ACGU':
        al = 'ACGU'
    else:
        al = 'ACDEFGHIKLMNPQRSTVWY'
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
    # Depending on the version of MEME, this line either like like
    #    MOTIF  1        width =  19  sites =   3  llr = 43  E-value = 6.9e-002
    # or like
    #    MOTIF  1 MEME    width =  19  sites =   3  llr = 43  E-value = 6.9e-002
    # or in v 4.11.4 onwards
    #    MOTIF ATTATAAAAAAA MEME-1	width =  12  sites =   5  llr = 43  E-value = 1.9e-003
    words = line.split()
    assert words[0] == 'MOTIF'
    if words[2][:5] == 'MEME-':
        motif_number = int(words[2].split('-')[1])
    else:
        motif_number = int(words[1])
    if words[2].startswith('MEME'):
        key_values = words[3:]
    else:
        key_values = words[2:]
    keys = key_values[::3]
    equal_signs = key_values[1::3]
    values = key_values[2::3]
    assert keys == ['width', 'sites', 'llr', 'E-value']
    for equal_sign in equal_signs:
        assert equal_sign == '='
    length = int(values[0])
    num_occurrences = int(values[1])
    evalue = float(values[3])
    return motif_number, length, num_occurrences, evalue


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
        assert len(sequence) == length
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
