# Copyright 2011, 2012 by Andrew Sczesnak.  All rights reserved.
# Revisions Copyright 2011, 2017 by Peter Cock.  All rights reserved.
# Revisions Copyright 2014, 2015 by Adam Novak.  All rights reserved.
# Revisions Copyright 2015, 2017 by Blaise Li.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for the "maf" multiple alignment format.

The Multiple Alignment Format, described by UCSC, stores a series of
multiple alignments in a single file. It is suitable for whole-genome
to whole-genome alignments, metadata such as source chromosome, start
position, size, and strand can be stored.

See http://genome.ucsc.edu/FAQ/FAQformat.html#format5

You are expected to use this module via the Bio.Align functions.

Coordinates in the MAF format are defined in terms of zero-based start
positions (like Python) and aligning region sizes.

A minimal aligned region of length one and starting at first position in the
source sequence would have ``start == 0`` and ``size == 1``.

As we can see on this example, ``start + size`` will give one more than the
zero-based end position. We can therefore manipulate ``start`` and
``start + size`` as python list slice boundaries.
"""
import shlex
import itertools


from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord


class AlignmentWriter(interfaces.AlignmentWriter):
    """Accepts Alignment objects, writes a MAF file."""

    fmt = "MAF"

    def _write_trackline(self, stream, metadata):
        stream.write("track")
        for key, value in metadata.items():
            if key in ("name", "description", "frames"):
                pass
            elif key == "mafDot":
                if value not in ("on", "off"):
                    raise ValueError(
                        "mafDot value must be 'on' or 'off' (received '%s')" % value
                    )
            elif key == "visibility":
                if value not in ("dense", "pack", "full"):
                    raise ValueError(
                        "visibility value must be 'dense', 'pack', or 'full' (received '%s')"
                        % value
                    )
            elif key == "speciesOrder":
                value = " ".join(value)
            else:
                continue
            if " " in value:
                value = '"%s"' % value
            stream.write(f" {key}={value}")
        stream.write("\n")

    def write_header(self, stream, alignments):
        """Write the MAF header."""
        try:
            metadata = alignments.metadata
        except AttributeError:
            metadata = {"MAF Version": "1"}
        track_keys = (
            "name",
            "description",
            "frames",
            "mafDot",
            "visibility",
            "speciesOrder",
        )
        for key in track_keys:
            if key in metadata:
                self._write_trackline(stream, metadata)
                break
        stream.write("##maf")
        for key, value in metadata.items():
            if key in track_keys:
                continue
            if key == "Comments":
                continue
            if key == "MAF Version":
                if value != "1":
                    raise ValueError("MAF version must be 1")
                key = "version"
            elif key == "Scoring":
                key = "scoring"
            elif key == "Program":
                key = "program"
            else:
                raise ValueError("Unexpected key '%s' for header" % key)
            stream.write(f" {key}={value}")
        stream.write("\n")
        comments = metadata.get("Comments")
        if comments is not None:
            for comment in comments:
                stream.write(f"# {comment}\n")
            stream.write("\n")

    def _format_score_line(self, alignment, annotations):
        try:
            score = alignment.score
        except AttributeError:
            line = "a"
        else:
            line = f"a score={score:.6f}"
        value = annotations.get("pass")
        if value is not None:
            line += f" pass={value}"
        return line + "\n"

    def format_alignment(self, alignment):
        """Return a string with a single alignment formatted as a MAF block."""
        if not isinstance(alignment, Alignment):
            raise TypeError("Expected an Alignment object")
        try:
            alignment_annotations = alignment.annotations
        except AttributeError:
            alignment_annotations = {}
        lines = []
        line = self._format_score_line(alignment, alignment_annotations)
        lines.append(line)
        name_width = 0
        start_width = 0
        size_width = 0
        length_width = 0
        n = len(alignment.sequences)
        for i in range(n):
            record = alignment.sequences[i]
            coordinates = alignment.coordinates[i]
            try:
                name = record.id
            except AttributeError:
                name = "sequence_%d" % i
            start = coordinates[0]
            end = coordinates[-1]
            length = len(record)
            if start < end:
                size = end - start
            else:
                size = start - end
                start = length - start
            name_width = max(name_width, len(name))
            start_width = max(start_width, len(str(start)))
            size_width = max(size_width, len(str(size)))
            length_width = max(length_width, len(str(length)))
        for i, empty in enumerate(alignment_annotations.get("empty", [])):
            record, segment, status = empty
            try:
                name = record.id
            except AttributeError:
                name = "sequence_%d" % (i + n)
            name_width = max(name_width, len(name))
            start, end = segment
            length = len(record.seq)
            if start <= end:
                size = end - start
            else:
                size = start - end
                start = length - start
            start_width = max(start_width, len(str(start)))
            size_width = max(size_width, len(str(size)))
            length_width = max(length_width, len(str(length)))
        quality_width = name_width + start_width + size_width + length_width + 5
        for i in range(n):
            record = alignment.sequences[i]
            coordinates = alignment.coordinates[i]
            try:
                record_id = record.id
            except AttributeError:
                record_id = "sequence_%d" % i
            start = coordinates[0]
            end = coordinates[-1]
            length = len(record)
            if start < end:
                size = end - start
                strand = "+"
            else:
                size = start - end
                start = length - start
                strand = "-"
            text = alignment[i]
            name = record_id.ljust(name_width)
            start = str(start).rjust(start_width)
            size = str(size).rjust(size_width)
            length = str(length).rjust(length_width)
            line = f"s {name} {start} {size} {strand} {length} {text}\n"
            lines.append(line)
            try:
                annotations = record.annotations
            except AttributeError:
                annotations = None
            if annotations is not None:
                quality = annotations.get("quality")
                if quality is not None:
                    gapped_quality = ""
                    j = 0
                    for letter in text:
                        if letter == "-":
                            gapped_quality += "-"
                        else:
                            gapped_quality += quality[j]
                            j += 1
                    name = record_id.ljust(quality_width)
                    line = f"q {name} {gapped_quality}\n"
                    lines.append(line)
                try:
                    leftStatus = annotations["leftStatus"]
                    leftCount = annotations["leftCount"]
                    rightStatus = annotations["rightStatus"]
                    rightCount = annotations["rightCount"]
                except KeyError:
                    pass
                else:
                    name = record_id.ljust(name_width)
                    line = f"i {name} {leftStatus} {leftCount} {rightStatus} {rightCount}\n"
                    lines.append(line)
        for i, empty in enumerate(alignment_annotations.get("empty", [])):
            record, segment, status = empty
            try:
                name = record.id
            except AttributeError:
                name = "sequence_%d" % (i + n)
            name = name.ljust(name_width)
            start, end = segment
            length = len(record.seq)
            if start <= end:
                size = end - start
                strand = "+"
            else:
                size = start - end
                start = length - start
                strand = "-"
            start = str(start).rjust(start_width)
            size = str(size).rjust(size_width)
            length = str(length).rjust(length_width)
            line = f"e {name} {start} {size} {strand} {length} {status}\n"
            lines.append(line)
        lines.append("\n")
        return "".join(lines)


class AlignmentIterator(interfaces.AlignmentIterator):
    """Alignment iterator for Multiple Alignment Format files.

    The file may contain multiple concatenated alignments, which are loaded
    and returned incrementally.

    File meta-data are stored in the ``.metadata`` attribute of the returned
    iterator.  Alignment annotations are stored in the ``.annotations``
    attribute of the ``Alignment`` object, except for the alignment score,
    which is stored as an attribute.  Sequence information of empty parts in
    the alignment block (sequences that connect the previous alignment block to
    the next alignment block, but do not align to the current alignment block)
    is stored in the alignment annotations under the ``"empty"`` key.
    Annotations specific to each line in the alignment are stored in the
    ``.annotations`` attribute of the corresponding sequence record.
    """

    fmt = "MAF"

    status_characters = ("C", "I", "N", "n", "M", "T")
    empty_status_characters = ("C", "I", "M", "n")

    def _read_header(self, stream):
        metadata = {}
        line = next(stream)
        if line.startswith("track "):
            words = shlex.split(line)
            for word in words[1:]:
                key, value = word.split("=")
                if key in ("name", "description", "frames"):
                    pass
                elif key == "mafDot":
                    if value not in ("on", "off"):
                        raise ValueError(
                            "Variable mafDot in track line has unexpected value '%s'"
                            % value
                        )
                elif key == "visibility":
                    if value not in ("dense", "pack", "full"):
                        raise ValueError(
                            "Variable visibility in track line has unexpected value '%s'"
                            % value
                        )
                elif key == "speciesOrder":
                    value = value.split()
                else:
                    raise ValueError("Unexpected variable '%s' in track line" % key)
                metadata[key] = value
            line = next(stream)
        words = line.split()
        if words[0] != "##maf":
            raise ValueError("header line does not start with ##maf")
        for word in words[1:]:
            key, value = word.split("=")
            if key == "version":
                key = "MAF Version"
            elif key == "scoring":
                key = "Scoring"
            elif key == "program":
                key = "Program"
            else:
                raise ValueError("Unexpected variable '%s' in header line" % key)
            metadata[key] = value
        if metadata.get("MAF Version") != "1":
            raise ValueError("MAF version must be 1")
        comments = []
        for line in stream:
            if line.strip():
                if not line.startswith("#"):
                    assert line.startswith("a")
                    self._aline = line
                    break
                comment = line[1:].strip()
                comments.append(comment)
        else:
            self._close()
        if comments:
            metadata["Comments"] = comments
        self.metadata = metadata

    def _read_next_alignment(self, stream):
        aline = self._aline
        if aline is None:
            return
        alignment = self._create_alignment(aline, stream)
        return alignment

    def _create_alignment(self, aline, stream):
        records = []
        strands = []
        aligned_sequences = []
        annotations = {}
        words = aline[1:].split()
        for word in words:
            key, value = word.split("=")
            if key == "score":
                score = float(value)
            elif key == "pass":
                value = int(value)
                if value <= 0:
                    raise ValueError("pass value must be positive (found %d)" % value)
                annotations["pass"] = value
            else:
                raise ValueError("Unknown annotation variable '%s'" % key)

        for line in stream:
            if line.startswith("#"):
                continue
            elif line.startswith("a"):
                self._aline = line
                break
            elif line.startswith("s "):
                words = line.strip().split()
                if len(words) != 7:
                    raise ValueError(
                        "Error parsing alignment - 's' line must have 7 fields"
                    )
                src = words[1]
                start = int(words[2])
                size = int(words[3])
                strand = words[4]
                srcSize = int(words[5])
                text = words[6]
                for gap_char in ".=_":
                    text = text.replace(gap_char, "-")
                aligned_sequences.append(text)
                sequence = text.replace("-", "")
                if len(sequence) != size:
                    raise ValueError(
                        "sequence size is incorrect (found %d, expected %d)"
                        % (len(sequence), size)
                    )
                if strand == "-":
                    sequence = reverse_complement(sequence)
                    start = srcSize - start - size
                seq = Seq({start: sequence}, length=srcSize)
                record = SeqRecord(seq, id=src, name="", description="")
                records.append(record)
                strands.append(strand)
            elif line.startswith("i "):
                words = line.strip().split()
                assert len(words) == 6
                assert words[1] == src  # from the previous "s" line
                leftStatus = words[2]
                leftCount = int(words[3])
                rightStatus = words[4]
                rightCount = int(words[5])
                assert leftStatus in AlignmentIterator.status_characters
                assert rightStatus in AlignmentIterator.status_characters
                record.annotations["leftStatus"] = leftStatus
                record.annotations["leftCount"] = leftCount
                record.annotations["rightStatus"] = rightStatus
                record.annotations["rightCount"] = rightCount
            elif line.startswith("e"):
                words = line[1:].split()
                assert len(words) == 6
                src = words[0]
                start = int(words[1])
                size = int(words[2])
                strand = words[3]
                srcSize = int(words[4])
                status = words[5]
                assert status in AlignmentIterator.empty_status_characters
                sequence = Seq(None, length=srcSize)
                record = SeqRecord(sequence, id=src, name="", description="")
                end = start + size
                if strand == "+":
                    segment = (start, end)
                else:
                    segment = (srcSize - start, srcSize - end)
                empty = (record, segment, status)
                annotation = annotations.get("empty")
                if annotation is None:
                    annotation = []
                    annotations["empty"] = annotation
                annotation.append(empty)
            elif line.startswith("q "):
                words = line.strip().split()
                assert len(words) == 3
                assert words[1] == src  # from the previous "s" line
                value = words[2].replace("-", "")
                record.annotations["quality"] = value
            elif not line.strip():
                # reached the end of the alignment, but keep reading until we
                # find the next alignment
                continue
            else:
                raise ValueError(f"Error parsing alignment - unexpected line:\n{line}")
        else:
            self._aline = None
        coordinates = Alignment.infer_coordinates(aligned_sequences)
        for record, strand, row in zip(records, strands, coordinates):
            if strand == "-":
                row[:] = row[-1] - row[0] - row
            start = record.seq.defined_ranges[0][0]
            row += start
        alignment = Alignment(records, coordinates)
        if annotations is not None:
            alignment.annotations = annotations
        if score is not None:
            alignment.score = score
        return alignment
