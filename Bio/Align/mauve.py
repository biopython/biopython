# Copyright 2015-2015 by Eric Rasche.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for "xmfa" output from Mauve/ProgressiveMauve.

You are expected to use this module via the Bio.Align functions.
"""

from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord


class AlignmentWriter(interfaces.AlignmentWriter):
    """Mauve xmfa alignment writer."""

    fmt = "Mauve"

    def __init__(self, target, metadata=None, identifiers=None):
        """Create an AlignmentWriter object.

        Arguments:
         - target       - output stream or file name
         - metadata     - metadata to be included in the output. If metadata
                          is None, then the alignments object to be written
                          must have an attribute `metadata`.
         - identifiers  - list of the IDs of the sequences included in the
                          alignment. Sequences will be numbered according to
                          their index in this list. If identifiers is None,
                          then the alignments object to be written must have
                          an attribute `identifiers`.
        """
        super().__init__(target)
        self._metadata = metadata
        self._identifiers = identifiers

    def write_header(self, alignments):
        """Write the file header to the output file."""
        stream = self.stream
        metadata = self._metadata
        format_version = metadata.get("FormatVersion", "Mauve1")
        line = f"#FormatVersion {format_version}\n"
        stream.write(line)
        identifiers = self._identifiers
        filename = metadata.get("File")
        if filename is None:
            # sequences came from separate files
            for index, filename in enumerate(identifiers):
                number = index + 1
                line = f"#Sequence{number}File\t{filename}\n"
                stream.write(line)
                line = f"#Sequence{number}Format\tFastA\n"
                stream.write(line)
        else:
            # sequences came from one combined file
            for number, identifier in enumerate(identifiers):
                assert number == int(identifier)
                number += 1
                line = f"#Sequence{number}File\t{filename}\n"
                stream.write(line)
                line = f"#Sequence{number}Entry\t{number}\n"
                stream.write(line)
                line = f"#Sequence{number}Format\tFastA\n"
                stream.write(line)
        backbone_file = metadata.get("BackboneFile")
        if backbone_file is not None:
            line = f"#BackboneFile\t{backbone_file}\n"
            stream.write(line)

    def write_file(self, alignments):
        """Write a file with the alignments, and return the number of alignments.

        alignments - A Bio.Align.mauve.AlignmentIterator object.
        """
        metadata = self._metadata
        if metadata is None:
            try:
                metadata = alignments.metadata
            except AttributeError:
                raise ValueError("alignments do not have an attribute `metadata`")
            else:
                self._metadata = metadata
        identifiers = self._identifiers
        if identifiers is None:
            try:
                identifiers = alignments.identifiers
            except AttributeError:
                raise ValueError("alignments do not have an attribute `identifiers`")
            else:
                self._identifiers = identifiers
        count = interfaces.AlignmentWriter.write_file(self, alignments)
        return count

    def format_alignment(self, alignment):
        """Return a string with a single alignment in the Mauve format."""
        metadata = self._metadata
        n, m = alignment.shape

        if n == 0:
            raise ValueError("Must have at least one sequence")
        if m == 0:
            raise ValueError("Non-empty sequences are required")

        filename = metadata.get("File")
        lines = []
        for i in range(n):
            identifier = alignment.sequences[i].id
            start = alignment.coordinates[i, 0]
            end = alignment.coordinates[i, -1]
            if start <= end:
                strand = "+"
            else:
                strand = "-"
                start, end = end, start
            if start == end:
                assert start == 0
            else:
                start += 1  # switch to 1-based counting
            sequence = alignment[i]
            if filename is None:
                number = (
                    self._identifiers.index(identifier) + 1
                )  # Switch to 1-based counting
                line = f"> {number}:{start}-{end} {strand} {identifier}\n"
            else:
                number = int(identifier) + 1  # Switch to 1-based counting
                line = f"> {number}:{start}-{end} {strand} {filename}\n"
            lines.append(line)
            line = f"{sequence}\n"
            lines.append(line)
        lines.append("=\n")
        return "".join(lines)


class AlignmentIterator(interfaces.AlignmentIterator):
    """Mauve xmfa alignment iterator."""

    fmt = "Mauve"

    def _read_header(self, stream):
        metadata = {}
        prefix = "Sequence"
        suffixes = ("File", "Entry", "Format")
        id_info = {}
        for suffix in suffixes:
            id_info[suffix] = []
        for line in stream:
            if not line.startswith("#"):
                self._line = line.strip()
                break
            key, value = line[1:].split()
            if key.startswith(prefix):
                for suffix in suffixes:
                    if key.endswith(suffix):
                        break
                else:
                    raise ValueError("Unexpected keyword '%s'" % key)
                if suffix == "Entry":
                    value = int(value) - 1  # Switch to 0-based counting
                seq_num = int(key[len(prefix) : -len(suffix)])
                id_info[suffix].append(value)
                assert seq_num == len(id_info[suffix])  # Mauve uses 1-based counting
            else:
                metadata[key] = value.strip()
        else:
            if not metadata:
                raise ValueError("Empty file.") from None
        if len(set(id_info["File"])) == 1:
            # A single file containing all sequences was provided as input;
            # store the file name once, and use the entry number as ID
            metadata["File"] = id_info["File"][0]
            self.identifiers = [str(entry) for entry in id_info["Entry"]]
        else:
            assert len(set(id_info["File"])) == len(id_info["File"])
            # Separate files for each of the sequences were provided as input;
            # use the sequence file as ID
            self.identifiers = id_info["File"]
        self.metadata = metadata

    def _parse_description(self, line):
        assert line.startswith(">")
        locus, strand, comments = line[1:].split(None, 2)
        seq_num, start_end = locus.split(":")
        seq_num = int(seq_num) - 1  # python counting
        identifier = self.identifiers[seq_num]
        assert strand in "+-"
        start, end = start_end.split("-")
        start = int(start)
        end = int(end)
        if start == 0:
            assert end == 0  # unaligned sequence
        else:
            start -= 1  # python counting
        return (identifier, start, end, strand, comments)

    def _read_next_alignment(self, stream):
        descriptions = []
        seqs = []

        try:
            line = self._line
        except AttributeError:
            pass
        else:
            del self._line
            description = self._parse_description(line)
            identifier, start, end, strand, comments = description
            descriptions.append(description)
            seqs.append("")

        for line in stream:
            line = line.strip()
            if line.startswith("="):
                # There may be more data, but we've reached the end of this
                # alignment
                coordinates = Alignment.infer_coordinates(seqs)
                records = []
                for index, (description, seq) in enumerate(zip(descriptions, seqs)):
                    identifier, start, end, strand, comments = description
                    seq = seq.replace("-", "")
                    assert len(seq) == end - start
                    if strand == "+":
                        pass
                    elif strand == "-":
                        seq = reverse_complement(seq, inplace=False)
                        coordinates[index, :] = len(seq) - coordinates[index, :]
                    else:
                        raise ValueError("Unexpected strand '%s'" % strand)
                    coordinates[index] += start
                    if start == 0:
                        seq = Seq(seq)
                    else:
                        seq = Seq({start: seq}, length=end)
                    record = SeqRecord(seq, id=identifier, description=comments)
                    records.append(record)

                return Alignment(records, coordinates)

            elif line.startswith(">"):
                description = self._parse_description(line)
                identifier, start, end, strand, comments = description
                descriptions.append(description)
                seqs.append("")
            else:
                seqs[-1] += line
