"""Bio.SeqIO support for the Graphical Fragment Assembly format."""

import warnings
import hashlib

from Bio import BiopythonWarning
from Bio.File import as_handle
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def _check_tags(seq, tags):
    """Check a segment line's tags for inconsistencies (PRIVATE)."""
    for tag in tags:
        if tag[:2] == "LN":
            # Sequence length
            if int(tag[5:]) != len(seq):
                warnings.warn(
                    f"Segment line has incorrect length. Expected {tag[5:]} but got {len(seq)}.",
                    BiopythonWarning,
                )
        if tag[:2] == "SH":
            # SHA256 checksum
            checksum = hashlib.sha256(str(seq).encode()).hexdigest()
            if checksum.upper() != tag[5:]:
                warnings.warn(
                    f"Segment line has incorrect checksum. Expected {tag[5:]} but got {checksum}.",
                    BiopythonWarning,
                )


def GfaIterator(source):
    """Parser for GFA 1.x files.

    Documentation: https://gfa-spec.github.io/GFA-spec/GFA1.html
    """
    with as_handle(source) as handle:
        for line in handle:
            if line == "\n":
                warnings.warn("GFA data has a blank line.", BiopythonWarning)

            fields = line.strip("\n").split("\t")
            if fields[0] != "S":
                continue
            if len(fields) < 3:
                raise ValueError(
                    f"Segment line must have name and sequence fields: {line}."
                )

            if fields[2] == "*":
                # no sequence data
                continue
            seq = Seq(fields[2])
            _check_tags(seq, fields[3:])

            yield SeqRecord(seq, name=fields[1])


def Gfa2Iterator(source):
    """Parser for GFA 2.0 files.

    Documentation for version 2: https://gfa-spec.github.io/GFA-spec/GFA2.html
    """
    with as_handle(source) as handle:
        for line in handle:
            if line == "\n":
                warnings.warn("GFA data has a blank line.", BiopythonWarning)

            fields = line.strip("\n").split("\t")
            if fields[0] != "S":
                continue
            if len(fields) < 4:
                raise ValueError(
                    f"Segment line must have name, length, and sequence fields: {line}."
                )
            try:
                int(fields[2])
            except ValueError:
                raise ValueError(
                    f"Segment line must have an integer length: {line}."
                ) from None

            if fields[3] == "*":
                # no sequence data
                continue
            seq = Seq(fields[3])
            _check_tags(seq, fields[3:])

            yield SeqRecord(seq, name=fields[1])
