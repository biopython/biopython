# Copyright 2018 by Ariel Aptekmann.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Module for the support of MEME minimal motif format."""

from Bio import motifs


def read(handle):
    """Parse the text output of the MEME program into a meme.Record object.

    Examples
    --------
    >>> from Bio.motifs import minimal
    >>> with open("motifs/meme.out") as f:
    ...     record = minimal.read(f)
    ...
    >>> for motif in record:
    ...     print(motif.name, motif.evalue)
    ...
    1 1.1e-22

    You can access individual motifs in the record by their index or find a motif
    by its name:

    >>> from Bio import motifs
    >>> with open("motifs/minimal_test.meme") as f:
    ...     record = motifs.parse(f, 'minimal')
    ...
    >>> motif = record[0]
    >>> print(motif.name)
    KRP
    >>> motif = record['IFXA']
    >>> print(motif.name)
    IFXA

    This function wont retrieve instances, as there are none in minimal meme format.

    """
    motif_number = 0
    record = Record()
    _read_version(record, handle)
    _read_alphabet(record, handle)
    _read_background(record, handle)

    while True:
        for line in handle:
            if line.startswith("MOTIF"):
                break
        else:
            return record
        name = line.split()[1]
        motif_number += 1
        length, num_occurrences, evalue = _read_motif_statistics(handle)
        counts = _read_lpm(handle, num_occurrences)
        # {'A': 0.25, 'C': 0.25, 'T': 0.25, 'G': 0.25}
        motif = motifs.Motif(alphabet=record.alphabet, counts=counts)
        motif.background = record.background
        motif.length = length
        motif.num_occurrences = num_occurrences
        motif.evalue = evalue
        motif.name = name
        record.append(motif)
        assert len(record) == motif_number
    return record


class Record(list):
    """Class for holding the results of a minimal MEME run."""

    def __init__(self):
        """Initialize record class values."""
        self.version = ""
        self.datafile = ""
        self.command = ""
        self.alphabet = None
        self.background = {}
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


def _read_background(record, handle):
    """Read background letter frequencies (PRIVATE)."""
    for line in handle:
        if line.startswith("Background letter frequencies"):
            break
    else:
        raise ValueError(
            "Improper input file. File should contain a line starting background frequencies."
        )
    try:
        line = next(handle)
    except StopIteration:
        raise ValueError(
            "Unexpected end of stream: Expected to find line starting background frequencies."
        )
    line = line.strip()
    ls = line.split()
    A, C, G, T = float(ls[1]), float(ls[3]), float(ls[5]), float(ls[7])
    record.background = {"A": A, "C": C, "G": G, "T": T}


def _read_version(record, handle):
    """Read MEME version (PRIVATE)."""
    for line in handle:
        if line.startswith("MEME version"):
            break
    else:
        raise ValueError(
            "Improper input file. File should contain a line starting MEME version."
        )
    line = line.strip()
    ls = line.split()
    record.version = ls[2]


def _read_alphabet(record, handle):
    """Read alphabet (PRIVATE)."""
    for line in handle:
        if line.startswith("ALPHABET"):
            break
    else:
        raise ValueError(
            "Unexpected end of stream: Expected to find line starting with 'ALPHABET'"
        )
    if not line.startswith("ALPHABET= "):
        raise ValueError("Line does not start with 'ALPHABET':\n%s" % line)
    line = line.strip().replace("ALPHABET= ", "")
    if line == "ACGT":
        al = "ACGT"
    else:
        al = "ACDEFGHIKLMNPQRSTVWY"
    record.alphabet = al


def _read_lpm(handle, num_occurrences):
    """Read letter probability matrix (PRIVATE)."""
    counts = [[], [], [], []]
    for line in handle:
        freqs = line.split()
        if len(freqs) != 4:
            break
        counts[0].append(round(float(freqs[0]) * num_occurrences))
        counts[1].append(round(float(freqs[1]) * num_occurrences))
        counts[2].append(round(float(freqs[2]) * num_occurrences))
        counts[3].append(round(float(freqs[3]) * num_occurrences))
    c = {}
    c["A"] = counts[0]
    c["C"] = counts[1]
    c["G"] = counts[2]
    c["T"] = counts[3]
    return c


def _read_motif_statistics(handle):
    """Read motif statistics (PRIVATE)."""
    # minimal :
    #      letter-probability matrix: alength= 4 w= 19 nsites= 17 E= 4.1e-009
    for line in handle:
        if line.startswith("letter-probability matrix:"):
            break
    num_occurrences = int(line.split("nsites=")[1].split()[0])
    length = int(line.split("w=")[1].split()[0])
    evalue = float(line.split("E=")[1].split()[0])
    return length, num_occurrences, evalue


def _read_motif_name(handle):
    """Read motif name (PRIVATE)."""
    for line in handle:
        if "sorted by position p-value" in line:
            break
    else:
        raise ValueError("Unexpected end of stream: Failed to find motif name")
    line = line.strip()
    words = line.split()
    name = " ".join(words[0:2])
    return name
