# Copyright 2021 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for the "psl" multiple alignment format.

The Pattern Space Layout (PSL foramt, described by UCSC, stores a series
of pairwise alignments in a single file. Typically they are used for
transcript to genome alignments. PSL files stored the alignment positions
and alignment scores, but do not store the aligned sequences.

See http://genome.ucsc.edu/FAQ/FAQformat.html#format5

You are expected to use this module via the Bio.Align functions.

Coordinates in the PSL format are defined in terms of zero-based start
positions (like Python) and aligning region sizes.

A minimal aligned region of length one and starting at first position in the
source sequence would have ``start == 0`` and ``size == 1``.

As we can see on this example, ``start + size`` will give one more than the
zero-based end position. We can therefore manipulate ``start`` and
``start + size`` as python list slice boundaries.
"""
from itertools import chain
import numpy


from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class AlignmentWriter(interfaces.AlignmentWriter):
    """Accepts Alignment objects, writes a MAF file."""

    def _write_trackline(self, metadata):
        stream = self.stream
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

    def write_header(self, alignments):
        """Write the MAF header."""
        stream = self.stream
        metadata = alignments.metadata
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
                self._write_trackline(metadata)
                break
        stream.write("##maf")
        for key, value in metadata.items():
            if key in track_keys:
                continue
            if key == "comments":
                continue
            if key not in ("version", "scoring", "program"):
                raise ValueError("Unexpected key '%s' for header" % key)
            if key == "version" and value != "1":
                raise ValueError("MAF version must be 1")
            stream.write(f" {key}={value}")
        stream.write("\n")
        comments = metadata.get("comments")
        if comments is not None:
            for comment in comments:
                stream.write(f"# {comment}\n")
            stream.write("\n")

    def _write_score_line(self, alignment, annotations):
        stream = self.stream
        stream.write("a")
        try:
            score = alignment.score
        except AttributeError:
            pass
        else:
            stream.write(f" score={score:.6f}")
        if annotations is not None:
            value = annotations.get("pass")
            if value is not None:
                stream.write(f" pass={value}")
        stream.write("\n")

    def write_alignment(self, alignment):
        """Write a complete alignment to a MAF block."""
        if not isinstance(alignment, Alignment):
            raise TypeError("Expected an Alignment object")
        try:
            alignment_annotations = alignment.annotations
        except AttributeError:
            alignment_annotations = None
        self._write_score_line(alignment, alignment_annotations)
        name_width = 0
        start_width = 0
        size_width = 0
        length_width = 0
        stream = self.stream
        n = len(alignment.sequences)
        for i in range(n):
            record = alignment.sequences[i]
            coordinates = alignment.coordinates[i]
            name = record.id
            start = coordinates[0]
            end = coordinates[-1]
            length = len(record.seq)
            if start < end:
                size = end - start
            else:
                size = start - end
                start = length - start
            name_width = max(name_width, len(name))
            start_width = max(start_width, len(str(start)))
            size_width = max(size_width, len(str(size)))
            length_width = max(length_width, len(str(length)))
        for empty in alignment_annotations.get("empty", []):
            record, segment, status = empty
            name = record.id
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
            name = record.id
            start = coordinates[0]
            end = coordinates[-1]
            length = len(record.seq)
            if start < end:
                size = end - start
                strand = "+"
            else:
                size = start - end
                start = length - start
                strand = "-"
            text = alignment[i]
            name = record.id.ljust(name_width)
            start = str(start).rjust(start_width)
            size = str(size).rjust(size_width)
            length = str(length).rjust(length_width)
            line = f"s {name} {start} {size} {strand} {length} {text}\n"
            stream.write(line)
            try:
                annotations = record.annotations
            except AttributeError:
                annotations = None
            if annotations is not None:
                quality = annotations.get("quality")
                if quality is not None:
                    gapped_quality = ""
                    i = 0
                    for letter in text:
                        if letter == "-":
                            gapped_quality += "-"
                        else:
                            gapped_quality += quality[i]
                            i += 1
                    name = record.id.ljust(quality_width)
                    line = f"q {name} {gapped_quality}\n"
                    stream.write(line)
                try:
                    leftStatus = annotations["leftStatus"]
                    leftCount = annotations["leftCount"]
                    rightStatus = annotations["rightStatus"]
                    rightCount = annotations["rightCount"]
                except KeyError:
                    pass
                else:
                    name = record.id.ljust(name_width)
                    line = f"i {name} {leftStatus} {leftCount} {rightStatus} {rightCount}\n"
                    stream.write(line)
        for empty in alignment_annotations.get("empty", []):
            record, segment, status = empty
            name = record.id.ljust(name_width)
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
            stream.write(line)
        stream.write("\n")


class AlignmentIterator(interfaces.AlignmentIterator):
    """Alignment iterator for Pattern Space Layout (PSL) files.

    Each line in the file contains one pairwise alignments, which are loaded
    and returned incrementally.  Alignment score information such as the number
    of matches and mismatches are stored as attributes of each alignment.
    """

    def __init__(self, source):
        """Create an AlignmentIterator object.

        Arguments:
         - source   - input data or file name

        """
        super().__init__(source, mode="t", fmt="PSL")
        stream = self.stream
        metadata = {}
        line = next(stream)
        if line.startswith("psLayout "):
            words = line.split()
            if words[1] != "version":
                raise ValueError("Unexpected word '%s' in header line" % words[1])
            self.metadata = {"version": words[2]}
            line = next(stream)
            line = next(stream)
            line = next(stream)
            line = next(stream)
            if line.lstrip("-").strip() != "":
                raise ValueError("End of header not found")
            self.line = None
        else:
            self.line = line

    @staticmethod
    def create_alignment(
        records,
        aligned_sequences,
        strands,
        annotations,
        column_annotations,
        score,
    ):
        """Create the Alignment object from the collected alignment data."""
        coordinates = Alignment.infer_coordinates(aligned_sequences)
        for record, strand, row in zip(records, strands, coordinates):
            if strand == "-":
                row[:] = row[-1] - row[0] - row
            start = record.seq.defined_ranges[0][0]
            row += start
        alignment = Alignment(records, coordinates)
        if annotations is not None:
            alignment.annotations = annotations
        if column_annotations is not None:
            alignment.column_annotations = column_annotations
        if score is not None:
            alignment.score = score
        return alignment

    def parse(self, stream):
        """Parse the next alignment from the stream."""
        if stream is None:
            raise StopIteration

        line = self.line
        self.line = None
        if line is not None:
            lines = chain([line], stream)
        else:
            lines = stream
        for line in lines:
            words = line.split()
            if len(words) != 21:
                raise ValueError("line has %d columns; expected 21" % len(words))
            strand = words[8]
            qName = words[9]
            qSize = int(words[10])
            tName = words[13]
            tSize = int(words[14])
            blockCount = int(words[17])
            blockSizes = [int(blockSize) for blockSize in words[18].rstrip(",").split(",")]
            qStarts = [int(start) for start in words[19].rstrip(",").split(",")]
            tStarts = [int(start) for start in words[20].rstrip(",").split(",")]
            if len(blockSizes) != blockCount:
                raise ValueError("Inconsistent number of blocks (%d found, expected %d)" % (len(blockSizes), blockCount))
            if len(qStarts) != blockCount:
                raise ValueError("Inconsistent number of query start positions (%d found, expected %d)" % (len(qStarts), blockCount))
            if len(tStarts) != blockCount:
                raise ValueError("Inconsistent number of target start positions (%d found, expected %d)" % (len(tStarts), blockCount))
            target_sequence = Seq(None, length=tSize)
            target_record = SeqRecord(target_sequence, id=tName)
            query_sequence = Seq(None, length=qSize)
            query_record = SeqRecord(query_sequence, id=qName)
            records = [target_record, query_record]
            blockSizes = numpy.array(blockSizes)
            qStarts = numpy.array(qStarts)
            tStarts = numpy.array(tStarts)
            if strand in ("++", "+-"):
                # protein sequence aligned against translated DNA sequence
                qStarts *= 3
                qSize *= 3
                blockSizes *= 3
            qPosition = qStarts[0]
            tPosition = tStarts[0]
            coordinates = [[tPosition, qPosition]]
            for blockSize, tStart, qStart in zip(blockSizes, tStarts, qStarts):
                if tStart != tPosition:
                    coordinates.append([tStart, qPosition])
                    tPosition = tStart
                if qStart != qPosition:
                    coordinates.append([tPosition, qStart])
                    qPosition = qStart
                tPosition += blockSize
                qPosition += blockSize
                coordinates.append([tPosition, qPosition])
            coordinates = numpy.array(coordinates).transpose()
            if strand == "-":
                coordinates[1, :] = qSize - coordinates[1, :]
            elif strand == "+-":
                coordinates[0, :] = tSize - coordinates[0, :]
            alignment = Alignment(records, coordinates)
            alignment.matches = int(words[0])
            alignment.misMatches = int(words[1])
            alignment.repMatches = int(words[2])
            alignment.nCount = int(words[3])
            qStart = int(words[11])
            qEnd = int(words[12])
            tStart = int(words[15])
            tEnd = int(words[16])
            if strand == "-":
                qStart, qEnd = qEnd, qStart
            if strand == "+-":
                tStart, tEnd = tEnd, tStart
            if strand in ("++", "+-"):
                # protein sequence aligned against translated DNA sequence
                qStart *= 3
                qEnd *= 3
                qSize *= 3
            if tStart != coordinates[0, 0]:
                raise ValueError("Inconsistent tStart found (%d, expected %d)" % (tStart, coordinates[0, 0]))
            if tEnd != coordinates[0, -1]:
                raise ValueError("Inconsistent tEnd found (%d, expected %d)" % (tEnd, coordinates[0, -1]))
            if qStart != coordinates[1, 0]:
                raise ValueError("Inconsistent qStart found (%d, expected %d)" % (qStart, coordinates[1, 0]))
            if qEnd != coordinates[1, -1]:
                raise ValueError("Inconsistent qEnd found (%d, expected %d)" % (qEnd, coordinates[1, -1]))
            if strand == "-":
                coordinates = coordinates.copy()
                coordinates[1, :] = qSize - coordinates[1, :]
            if strand == "+-":
                coordinates = coordinates.copy()
                coordinates[0, :] = tSize - coordinates[0, :]
            qNumInsert = 0
            qBaseInsert = 0
            tNumInsert = 0
            tBaseInsert = 0
            tStart, qStart = coordinates[:, 0]
            for tEnd, qEnd in coordinates[:, 1:].transpose():
                tCount = tEnd - tStart
                qCount = qEnd - qStart
                if tCount == 0:
                    if qStart > 0 and qEnd < qSize:
                        qNumInsert += 1
                        qBaseInsert += qCount
                    qStart = qEnd
                elif qCount == 0:
                    if tStart > 0 and tEnd < tSize:
                        tNumInsert += 1
                        tBaseInsert += tCount
                    tStart = tEnd
                else:
                    tStart = tEnd
                    qStart = qEnd
            if strand in ("++", "+-"):
                assert qBaseInsert % 3 == 0
                qBaseInsert //= 3
            if qNumInsert != int(words[4]):
                raise ValueError("Inconsistent qNumInsert found (%s, expected %d)" % (words[4], qNumInsert))
            if qBaseInsert != int(words[5]):
                raise ValueError("Inconsistent qBaseInsert found (%s, expected %d)" % (words[5], qBaseInsert))
            if tNumInsert != int(words[6]):
                raise ValueError("Inconsistent tNumInsert found (%s, expected %d)" % (words[6], tNumInsert))
            if tBaseInsert != int(words[7]):
                raise ValueError("Inconsistent tBaseInsert found (%s, expected %d)" % (words[7], tBaseInsert))
            yield alignment
