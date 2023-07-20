# Copyright 2003 by Bartek Wilczynski.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Parsing AlignACE output files."""

from Bio.motifs import Motif
from Bio.Align import Alignment
from Bio.Seq import Seq


class Record(list):
    """AlignACE record (subclass of Python list)."""

    def __init__(self):
        """Initialize the class."""
        self.parameters = None


def read(handle):
    """Parse an AlignACE format handle as a Record object."""
    record = Record()
    line = next(handle)
    record.version = line.strip()
    line = next(handle)
    record.command = line.strip()
    mask = None
    number = None
    for line in handle:
        line = line.strip()
        if line == "":
            pass
        elif line[:4] == "Para":
            record.parameters = {}
        elif line[0] == "#":
            seq_name = line.split("\t")[1]
            record.sequences.append(seq_name)
        elif "=" in line:
            par_name, par_value = line.split("=")
            par_name = par_name.strip()
            par_value = par_value.strip()
            record.parameters[par_name] = par_value
        elif line[:5] == "Input":
            record.sequences = []
        elif line[:5] == "Motif":
            words = line.split()
            assert words[0] == "Motif"
            number = int(words[1])
            instances = []
        elif line[:3] == "MAP":
            alphabet = "ACGT"
            alignment = Alignment(instances)
            motif = Motif(alphabet, alignment)
            motif.score = float(line.split()[-1])
            motif.number = number
            motif.mask = mask
            record.append(motif)
        elif len(line.split("\t")) == 4:
            seq = Seq(line.split("\t")[0])
            instances.append(seq)
        elif "*" in line:
            mask = line.strip("\r\n")
        else:
            raise ValueError(line)
    return record
