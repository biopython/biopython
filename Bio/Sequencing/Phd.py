# Copyright 2004 by Cymon J. Cox and Frank Kauff.  All rights reserved.
# Copyright 2008 by Michiel de Hoon.  All rights reserved.
# Revisions copyright 2009 by Cymon J. Cox.  All rights reserved.
# Revisions copyright 2009 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Parser for PHD files output by PHRED and used by PHRAP and CONSED.

This module can be used directly which will return Record objects
which should contain all the original data in the file.

Alternatively, using Bio.SeqIO with the "phd" format will call this module
internally.  This will give SeqRecord objects for each contig sequence.
"""

from Bio import Seq
from Bio.Alphabet import generic_dna


CKEYWORDS = [
    "CHROMAT_FILE",
    "ABI_THUMBPRINT",
    "PHRED_VERSION",
    "CALL_METHOD",
    "QUALITY_LEVELS",
    "TIME",
    "TRACE_ARRAY_MIN_INDEX",
    "TRACE_ARRAY_MAX_INDEX",
    "TRIM",
    "TRACE_PEAK_AREA_RATIO",
    "CHEM",
    "DYE",
]


class Record(object):
    """Hold information from a PHD file."""

    def __init__(self):
        """Initialize the class."""
        self.file_name = ""
        self.comments = {}
        for kw in CKEYWORDS:
            self.comments[kw.lower()] = None
        self.sites = []
        self.seq = ""
        self.seq_trimmed = ""


def read(handle):
    """Read the next PHD record from the file, return it as a Record object.

    This function reads PHD file data line by line from the handle,
    and returns a single Record object.
    """
    for line in handle:
        if line.startswith("BEGIN_SEQUENCE"):
            record = Record()
            record.file_name = line[15:].rstrip()
            break
    else:
        return  # No record found

    for line in handle:
        if line.startswith("BEGIN_COMMENT"):
            break
    else:
        raise ValueError("Failed to find BEGIN_COMMENT line")

    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line == "END_COMMENT":
            break
        keyword, value = line.split(":", 1)
        keyword = keyword.lower()
        value = value.strip()
        if keyword in (
            "chromat_file",
            "phred_version",
            "call_method",
            "chem",
            "dye",
            "time",
            "basecaller_version",
            "trace_processor_version",
        ):
            record.comments[keyword] = value
        elif keyword in (
            "abi_thumbprint",
            "quality_levels",
            "trace_array_min_index",
            "trace_array_max_index",
        ):
            record.comments[keyword] = int(value)
        elif keyword == "trace_peak_area_ratio":
            record.comments[keyword] = float(value)
        elif keyword == "trim":
            first, last, prob = value.split()
            record.comments[keyword] = (int(first), int(last), float(prob))
    else:
        raise ValueError("Failed to find END_COMMENT line")

    for line in handle:
        if line.startswith("BEGIN_DNA"):
            break
    else:
        raise ValueError("Failed to find BEGIN_DNA line")

    for line in handle:
        if line.startswith("END_DNA"):
            break
        else:
            # Line is: "site quality peak_location"
            # Peak location is optional according to
            # David Gordon (the Consed author)
            parts = line.split()
            if len(parts) in [2, 3]:
                record.sites.append(tuple(parts))
            else:
                raise ValueError(
                    "DNA line must contain a base and quality "
                    "score, and optionally a peak location."
                )

    for line in handle:
        if line.startswith("END_SEQUENCE"):
            break
    else:
        raise ValueError("Failed to find END_SEQUENCE line")

    record.seq = Seq.Seq("".join(n[0] for n in record.sites), generic_dna)
    if record.comments["trim"] is not None:
        first, last = record.comments["trim"][:2]
        record.seq_trimmed = record.seq[first:last]

    return record


def parse(handle):
    """Iterate over a file yielding multiple PHD records.

    The data is read line by line from the handle. The handle can be a list
    of lines, an open file, or similar; the only requirement is that we can
    iterate over the handle to retrieve lines from it.

    Typical usage::

        records = parse(handle)
        for record in records:
            # do something with the record object

    """
    while True:
        record = read(handle)
        if not record:
            return
        yield record
