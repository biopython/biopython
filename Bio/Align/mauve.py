# Copyright 2015-2015 by Eric Rasche.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for "xmfa" output from Mauve/ProgressiveMauve.

You are expected to use this module via the Bio.Align functions.
"""

from Bio.Align import interfaces, Alignment
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord


class AlignmentWriter(interfaces.AlignmentWriter):
    """Mauve/XMFA alignment writer."""

    def write_header(self, alignments):
        """Write the file header to the output file."""
        stream = self.stream
        metadata = alignments.metadata
        format_version = metadata.get("FormatVersion", "Mauve1")
        line = f"#FormatVersion {format_version}\n"
        stream.write(line)
        alignment = alignments[0]
        if len(alignment) > 1:
            # this is a real alignment and includes all sequences
            names = [sequence.id for sequence in alignment.sequences]
        else:
            # this is a single sequence; file contains no real alignments
            # but lists the individual sequences only
            names = [alignment.sequences[0].id for alignment in alignments]
        name = names[0]
        try:
            filename, index = name.rsplit(":", 1)
        except ValueError:
            # sequences came from separate files
            for index, name in enumerate(names):
                number = index + 1
                line = f"#Sequence{number}File\t{name}\n"
                stream.write(line)
                line = f"#Sequence{number}Format\tFastA\n"
                stream.write(line)
            self.names = names
        else:
            # sequences came from one combined file
            for number in range(1, len(names) + 1):
                line = f"#Sequence{number}File\t{filename}\n"
                stream.write(line)
                line = f"#Sequence{number}Entry\t{number}\n"
                stream.write(line)
                line = f"#Sequence{number}Format\tFastA\n"
                stream.write(line)
        backbone_file = metadata.get("BackboneFile", None)
        if backbone_file is not None:
            line = f"#BackboneFile\t{backbone_file}\n"
            stream.write(line)

    def write_file(self, alignments):
        """Write a file with the alignments, and return the number of alignments.

        alignments - A Bio.Align.mauve.AlignmentIterator object.
        """

        class ListWithAttributes(list):
            pass

        try:
            metadata = alignments.metadata
        except AttributeError:
            metadata = {}
        alignments = ListWithAttributes(alignments)
        alignments.metadata = metadata
        count = interfaces.AlignmentWriter.write_file(self, alignments)
        return count

    def write_alignment(self, alignment):
        """Use this to write (another) single alignment to an open file."""
        n, m = alignment.shape

        if n == 0:
            raise ValueError("Must have at least one sequence")
        if m == 0:
            raise ValueError("Non-empty sequences are required")

        stream = self.stream
        for i in range(n):
            filename = alignment.sequences[i].id
            try:
                filename, number = filename.rsplit(":", 1)
            except ValueError:
                number = self.names.index(filename)
            else:
                number = int(number)
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
            number += 1  # switch to 1-based counting
            sequence = alignment[i]
            line = f"> {number}:{start}-{end} {strand} {filename}\n"
            stream.write(line)
            line = f"{sequence}\n"
            stream.write(line)
        stream.write("=\n")


class AlignmentIterator(interfaces.AlignmentIterator):
    """Mauve xmfa alignment iterator."""

    def __init__(self, source):
        """Create an AlignmentIterator object.

        Arguments:
         - source   - input data or file name

        """
        super().__init__(source, mode="t", fmt="Mauve")
        stream = self.stream
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
                seq_num = int(key[len(prefix) : -len(suffix)])
                id_info[suffix].append(value)
                assert seq_num == len(id_info[suffix])  # Mauve uses 1-based counting
            else:
                metadata[key] = value.strip()
        else:
            if not metadata:
                raise ValueError("Empty file.") from None
        if len(set(id_info["File"])) == 1:
            self._identifiers = []
            # Use filename + entry number as ID
            for filename, entry in zip(id_info["File"], id_info["Entry"]):
                entry = int(entry) - 1  # Use 0-based counting
                identifier = "%s:%d" % (filename, entry)
                self._identifiers.append(identifier)
        else:
            assert len(set(id_info["File"])) == len(id_info["File"])
            # All filenames are unique
            self._identifiers = id_info["File"]
        self.metadata = metadata

    def _parse_description(self, line):
        assert line.startswith(">")
        locus, strand, comments = line[1:].split(None, 2)
        seq_num, start_end = locus.split(":")
        seq_num = int(seq_num) - 1  # python counting
        identifier = self._identifiers[seq_num]
        assert strand in "+-"
        start, end = start_end.split("-")
        start = int(start)
        end = int(end)
        if start == 0:
            assert end == 0  # unaligned sequence
        else:
            start -= 1  # python counting
        return (identifier, start, end, strand, comments)

    def parse(self, stream):
        """Parse the next alignment from the stream."""
        if stream is None:
            return

        descriptions = []
        seqs = []

        line = self._line
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
                    length = end - start
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

                yield Alignment(records, coordinates)

                descriptions = []
                seqs = []
            elif line.startswith(">"):
                description = self._parse_description(line)
                identifier, start, end, strand, comments = description
                descriptions.append(description)
                seqs.append("")
            else:
                seqs[-1] += line
