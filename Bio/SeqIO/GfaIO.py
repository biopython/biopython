"""Bio.SeqIO support for the Graphical Fragment Assembly format.

This format is output by many assemblers and includes linkage information for
how the different sequences fit together, however, we just care about the
segment (sequence) information.

Documentation:
- Version 1.x: https://gfa-spec.github.io/GFA-spec/GFA1.html
- Version 2.0: https://gfa-spec.github.io/GFA-spec/GFA2.html
"""

import hashlib
import re
import warnings

from Bio import BiopythonWarning
from Bio.Seq import _UndefinedSequenceData
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


from .Interfaces import _TextIOSource
from .Interfaces import SequenceIterator


def _check_tags(seq, tags):
    """Check a segment line's tags for inconsistencies (PRIVATE)."""
    for tag in tags:
        if tag[:2] == "LN":
            # Sequence length
            if len(seq) == 0:
                # No sequence data, set the sequence length
                seq._data = _UndefinedSequenceData(int(tag[5:]))
            elif int(tag[5:]) != len(seq):
                warnings.warn(
                    f"Segment line has incorrect length. Expected {tag[5:]} but got {len(seq)}.",
                    BiopythonWarning,
                )
        elif tag[:2] == "SH":
            # SHA256 checksum
            checksum = hashlib.sha256(str(seq).encode()).hexdigest()
            if checksum.upper() != tag[5:]:
                warnings.warn(
                    f"Segment line has incorrect checksum. Expected {tag[5:]} but got {checksum}.",
                    BiopythonWarning,
                )


def _tags_to_annotations(tags):
    """Build an annotations dictionary from a list of tags (PRIVATE)."""
    annotations = {}
    for tag in tags:
        parts = tag.split(":")
        if len(parts) < 3:
            raise ValueError(f"Segment line has invalid tag: {tag}.")
        if re.fullmatch(r"[A-Za-z][A-Za-z0-9]", parts[0]) is None:
            warnings.warn(
                f"Tag has invalid name: {parts[0]}. Are they tab delimited?",
                BiopythonWarning,
            )
        parts[2] = ":".join(parts[2:])  # tag value may contain : characters
        annotations[parts[0]] = (parts[1], parts[2])

        # Check type of the tag and raise warning on a mismatch. These RegExs
        # are part of the 1.0 standard.
        if parts[1] not in "AifZJHB":
            warnings.warn(f"Tag has invalid type: {parts[1]}", BiopythonWarning)
        elif parts[1] == "A" and re.fullmatch(r"[!-~]", parts[2]) is None:
            warnings.warn(
                f"Tag has incorrect type. Expected printable character, got {parts[2]}.",
                BiopythonWarning,
            )
        elif parts[1] == "i" and re.fullmatch(r"[-+]?[0-9]+", parts[2]) is None:
            warnings.warn(
                f"Tag has incorrect type. Expected signed integer, got {parts[2]}.",
                BiopythonWarning,
            )
        elif (
            parts[1] == "f"
            and re.fullmatch(r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?", parts[2])
            is None
        ):
            warnings.warn(
                f"Tag has incorrect type. Expected float, got {parts[2]}.",
                BiopythonWarning,
            )
        elif parts[1] == "Z" and re.fullmatch(r"[ !-~]+", parts[2]) is None:
            warnings.warn(
                f"Tag has incorrect type. Expected printable string, got {parts[2]}.",
                BiopythonWarning,
            )
        elif parts[1] == "J" and re.fullmatch(r"[ !-~]+", parts[2]) is None:
            warnings.warn(
                f"Tag has incorrect type. Expected JSON excluding new-line and tab characters, got {parts[2]}.",
                BiopythonWarning,
            )
        elif parts[1] == "H" and re.fullmatch(r"[0-9A-F]+", parts[2]) is None:
            warnings.warn(
                f"Tag has incorrect type. Expected byte array in hex format, got {parts[2]}.",
                BiopythonWarning,
            )
        elif (
            parts[1] == "B"
            and re.fullmatch(
                r"[cCsSiIf](,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+", parts[2]
            )
            is None
        ):
            warnings.warn(
                f"Tag has incorrect type. Expected array of integers or floats, got {parts[2]}.",
                BiopythonWarning,
            )
    return annotations


class Gfa1Iterator(SequenceIterator):
    """Parser for GFA 1.x files.

    Documentation: https://gfa-spec.github.io/GFA-spec/GFA1.html
    """

    modes = "t"

    def __init__(
        self,
        source: _TextIOSource,
    ) -> None:
        """Iterate over a GFA file as SeqRecord objects.

        Arguments:
         - source - input stream opened in text mode, or a path to a file
        """
        super().__init__(source, fmt="GFA 1.0")

    def __next__(self):
        for line in self.stream:
            if line == "\n":
                warnings.warn("GFA data has a blank line.", BiopythonWarning)
                continue

            fields = line.strip("\n").split("\t")
            if fields[0] == "S":
                break
        else:
            raise StopIteration
        if len(fields) < 3:
            raise ValueError(
                f"Segment line must have name and sequence fields: {line}."
            )

        if fields[2] == "*":
            seq = Seq(None, length=0)
        else:
            seq = Seq(fields[2])

        tags = fields[3:]
        _check_tags(seq, tags)
        annotations = _tags_to_annotations(tags)

        return SeqRecord(seq, id=fields[1], name=fields[1], annotations=annotations)


class Gfa2Iterator(SequenceIterator):
    """Parser for GFA 2.0 files.

    Documentation for version 2: https://gfa-spec.github.io/GFA-spec/GFA2.html
    """

    modes = "t"

    def __init__(
        self,
        source: _TextIOSource,
    ) -> None:
        """Iterate over a GFA file as SeqRecord objects.

        Arguments:
         - source - input stream opened in text mode, or a path to a file
        """
        super().__init__(source, fmt="GFA 2.0")

    def __next__(self):
        for line in self.stream:
            if line == "\n":
                warnings.warn("GFA data has a blank line.", BiopythonWarning)
                continue

            fields = line.strip("\n").split("\t")
            if fields[0] == "S":
                break
        else:
            raise StopIteration
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
            seq = Seq(None, length=0)
        else:
            seq = Seq(fields[3])

        tags = fields[4:]
        _check_tags(seq, tags)
        annotations = _tags_to_annotations(tags)

        return SeqRecord(seq, id=fields[1], name=fields[1], annotations=annotations)
