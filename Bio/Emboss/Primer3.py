# Copyright 2008 Michiel de Hoon.
# Revisions copyright 2009 Leighton Pritchard.
# Revisions copyright 2010 Peter Cock.
# All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Code to parse output from the EMBOSS eprimer3 program.

As elsewhere in Biopython there are two input functions, read and parse,
for single record output and multi-record output. For primer3, a single
record object is created for each target sequence and may contain
multiple primers.

i.e. If you ran eprimer3 with a single target sequence, use the read
function. If you ran eprimer3 with multiple targets, use the parse
function to iterate over the retsults.
"""


# --- primer3


class Record:
    """Represent information from a primer3 run finding primers.

    Members:

        - primers  - list of Primer objects describing primer pairs for
          this target sequence.
        - comments - the comment line(s) for the record

    """

    def __init__(self):
        """Initialize the class."""
        self.comments = ""
        self.primers = []


class Primers:
    """A primer set designed by Primer3.

    Members:

        - size - length of product, note you can use len(primer) as an
          alternative to primer.size

        - forward_seq
        - forward_start
        - forward_length
        - forward_tm
        - forward_gc

        - reverse_seq
        - reverse_start
        - reverse_length
        - reverse_tm
        - reverse_gc

        - internal_seq
        - internal_start
        - internal_length
        - internal_tm
        - internal_gc

    """

    def __init__(self):
        """Initialize the class."""
        self.size = 0
        self.forward_seq = ""
        self.forward_start = 0
        self.forward_length = 0
        self.forward_tm = 0.0
        self.forward_gc = 0.0
        self.reverse_seq = ""
        self.reverse_start = 0
        self.reverse_length = 0
        self.reverse_tm = 0.0
        self.reverse_gc = 0.0
        self.internal_seq = ""
        self.internal_start = 0
        self.internal_length = 0
        self.internal_tm = 0.0
        self.internal_gc = 0.0

    def __len__(self):
        """Length of the primer product (i.e. product size)."""
        return self.size


def parse(handle):
    """Iterate over primer3 output as Bio.Emboss.Primer3.Record objects."""
    # Skip blank lines at head of file
    while True:
        line = handle.readline()
        if line.strip():
            break  # Starting a record

    # Read each record
    record = None
    primer = None
    while True:
        if line.startswith(("# EPRIMER3", "# PRIMER3")):
            # Record data
            if record is not None:
                yield record
            record = Record()
            record.comments += line
            primer = None
        elif line.startswith("#"):
            if (
                line.strip()
                != "#                      Start  Len   Tm     GC%   Sequence"
            ):
                record.comments += line
        elif not line.strip():
            pass
        elif line[5:19] == "PRODUCT SIZE: ":
            primer = Primers()
            primer.size = int(line[19:])
            record.primers.append(primer)
        elif line[5:19] == "FORWARD PRIMER":
            words = line.split()
            if not primer or primer.size == 0:
                primer = Primers()
                record.primers.append(primer)
            primer.forward_start = int(words[2])
            primer.forward_length = int(words[3])
            primer.forward_tm = float(words[4])
            primer.forward_gc = float(words[5])
            primer.forward_seq = words[6]
        elif line[5:19] == "REVERSE PRIMER":
            words = line.split()
            if not primer or primer.size == 0:
                primer = Primers()
                record.primers.append(primer)
            primer.reverse_start = int(words[2])
            primer.reverse_length = int(words[3])
            primer.reverse_tm = float(words[4])
            primer.reverse_gc = float(words[5])
            primer.reverse_seq = words[6]
        elif line[5:19] == "INTERNAL OLIGO":
            words = line.split()
            if not primer or primer.size == 0:
                primer = Primers()
                record.primers.append(primer)
            primer.internal_start = int(words[2])
            primer.internal_length = int(words[3])
            primer.internal_tm = float(words[4])
            primer.internal_gc = float(words[5])
            try:
                primer.internal_seq = words[6]
            except IndexError:  # eprimer3 reports oligo without sequence
                primer.internal_seq = ""
        try:
            line = next(handle)
        except StopIteration:
            break
    if record:
        yield record


def read(handle):
    """Parse primer3 output into a Bio.Emboss.Primer3.Record object.

    This is for when there is one and only one target sequence. If
    designing primers for multiple sequences, use the parse function.
    """
    iterator = parse(handle)
    try:
        record = next(iterator)
    except StopIteration:
        raise ValueError("No records found in handle") from None
    try:
        next(iterator)
        raise ValueError("More than one record found in handle")
    except StopIteration:
        pass
    return record
