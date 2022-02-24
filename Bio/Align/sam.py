# Copyright 2022 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for the "sam" pairwise alignment format.

The Sequence Alignment/Map (SAM) format, created by Heng Li and Richard Durbin
at the Wellcome Trust Sanger Institute, stores a series of alignments to the
genome in a single file. Typically they are used for next-generation sequencing
data. SAM files store the alignment positions for mapped sequences, and may
also store the aligned sequences and other information associated with the
sequence.

See http://www.htslib.org/ for more information.

You are expected to use this module via the Bio.Align functions.

Coordinates in the SAM format are defined in terms of one-based start
positions; the parser converts these to zero-based coordinates to be consistent
with Python and other alignment formats.
"""
from itertools import chain
import copy

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Please install numpy if you want to use Bio.Align. "
        "See http://www.numpy.org/"
    ) from None

from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq, reverse_complement, UndefinedSequenceError
from Bio.SeqRecord import SeqRecord


class Coordinates(numpy.ndarray):
    """numpy array subclass for alignment coordinates.

    This class uses a ``cigar`` attribute to store a BAM/SAM-style cigar
    representation of the alignment. The cigar can contain information about
    the alignment that is not represented in the coordinates. For example, the
    cigar string distinguishes between deletions from the reference and skipped
    regions (e.g. introns) from the reference.

    To ensure that the coordinates and the cigar remain consistent with each
    other, the coordinates are read-only by default, and the cigar is stored
    as an immutable tuple of tuples.
    """

    def __new__(cls, cigar, pos=0, strand="+"):
        """Create a new Coordinates instance."""
        query_pos = 0
        coordinates = [[pos, query_pos]]
        for operation, length in cigar:
            if operation in (0, 7, 8):  # M=X
                pos += length
                query_pos += length
            elif operation in (1, 4):  # IS
                query_pos += length
            elif operation in (2, 3):  # DN
                pos += length
            elif operation == 5:  # H
                # hard clipping (clipped sequences not present in sequence)
                continue
            elif operation == 6:  # P
                raise NotImplementedError("padding operator is not yet implemented")
            coordinates.append([pos, query_pos])
        shape = (len(coordinates), 2)
        obj = super().__new__(cls, shape, int)
        obj[:, :] = coordinates
        obj = obj.transpose()
        if strand == "-":
            obj[1, :] = query_pos - obj[1, :]
        obj.flags.writeable = False
        obj.cigar = tuple(cigar)
        return obj

    def copy(self):
        """Return a copy of the array as an unwriteable Coordinates object."""
        obj = super().copy()
        obj.flags.writeable = self.flags.writeable
        obj.cigar = self.cigar
        return obj


class AlignmentWriter(interfaces.AlignmentWriter):
    """Alignment file writer for the Sequence Alignment/Map (SAM) file format."""

    def __init__(self, target, md=False):
        """Create an AlignmentWriter object.

        Arguments:
         - md - If True, calculate the MD tag from the alignment and include it
                in the output.
                If False (default), do not include the MD tag in the output.

        """
        super().__init__(target, mode="w")
        self.md = md

    def write_header(self, alignments):
        """Write the SAM header."""
        try:
            metadata = alignments.metadata
        except AttributeError:
            metadata = {}
        try:
            targets = alignments.targets
        except AttributeError:
            targets = {}
        values = metadata.get("HD")
        if values is not None:
            # if HD is present, then VN is required and must come first
            fields = ["@HD", "VN:%s" % values["VN"]]
            for key, value in values.items():
                if key == "VN":
                    continue
                fields.append("%s:%s" % (key, value))
            line = "\t".join(fields) + "\n"
            self.stream.write(line)
        for rname, record in targets.items():
            assert rname == record.id
            fields = ["@SQ"]
            fields.append("SN:%s" % rname)
            length = len(record.seq)
            fields.append("LN:%d" % length)
            for key, value in record.annotations.items():
                if key == "alternate_locus":
                    fields.append("AH:%s" % value)
                elif key == "names":
                    fields.append("AN:%s" % ",".join(value))
                elif key == "assembly":
                    fields.append("AS:%s" % value)
                elif key == "MD5":
                    fields.append("M5:%s" % value)
                elif key == "species":
                    fields.append("SP:%s" % value)
                elif key == "topology":
                    assert value in ("linear", "circular")
                    fields.append("PP:%s" % value)
                elif key == "URI":
                    fields.append("UR:%s" % value)
                else:
                    fields.append("%s:%s" % (key[:2], value))
            try:
                description = record.description
            except AttributeError:
                pass
            else:
                if description != "<unknown description>":
                    fields.append("DS:%s" % description)
            line = "\t".join(fields) + "\n"
            self.stream.write(line)
        for tag, rows in metadata.items():
            if tag == "HD":  # already written
                continue
            for row in rows:
                fields = ["@" + tag]
                for key, value in row.items():
                    fields.append("%s:%s" % (key, value))
                line = "\t".join(fields) + "\n"
                self.stream.write(line)

    def format_alignment(self, alignment, md=None):
        """Return a string with a single alignment formatted as one SAM line."""
        if not isinstance(alignment, Alignment):
            raise TypeError("Expected an Alignment object")
        target, query = alignment.sequences
        try:
            qName = query.id
        except AttributeError:
            qName = "query"
            qual = "*"
        else:
            try:
                qual = query.letter_annotations["phred_quality"]
            except (AttributeError, KeyError):
                qual = "*"
            query = query.seq
        try:
            rName = target.id
        except AttributeError:
            rName = "target"
        else:
            target = target.seq
        coordinates = alignment.coordinates
        try:
            operations_lengths = coordinates.cigar
        except AttributeError:
            operations_lengths = None
        if coordinates[1, 0] < coordinates[1, -1]:  # mapped to forward strand
            flag = 0
        else:  # mapped to reverse strand
            flag = 16
            query = reverse_complement(query, inplace=False)
            coordinates = numpy.array(coordinates)
            coordinates[1, :] = len(query) - coordinates[1, :]
        try:
            query = bytes(query)
        except TypeError:  # string
            pass
        except UndefinedSequenceError:
            query = "*"
        else:
            query = str(query, "ASCII")
        if operations_lengths is None:
            # calculate the cigar from the alignment coordinates
            pos = None
            operations_lengths = []
            tSize = len(target)
            tStart, qStart = coordinates[:, 0]
            for tEnd, qEnd in coordinates[:, 1:].transpose():
                tCount = tEnd - tStart
                qCount = qEnd - qStart
                if tCount == 0:
                    length = qCount
                    if pos is None or tEnd == tSize:
                        operation = 4  # S; soft clipping
                    else:
                        operation = 1  # I; insertion to the reference
                    qStart = qEnd
                elif qCount == 0:
                    if tStart > 0 and tEnd < tSize:
                        length = tCount
                        operation = 2  # D; deletion from the reference
                    else:
                        operation = None
                    tStart = tEnd
                else:
                    if tCount != qCount:
                        raise ValueError("Unequal step sizes in alignment")
                    if pos is None:
                        pos = tStart
                    tStart = tEnd
                    qStart = qEnd
                    operation = 0  # M; alignment match
                    length = tCount
                if operation is not None:
                    operations_lengths.append([operation, length])
        else:
            # use the existing cigar
            tStart, qStart = coordinates[:, 0]
            for tEnd, qEnd in coordinates[:, 1:].transpose():
                if tStart < tEnd and qStart < qEnd:
                    pos = tStart
                    break
                tStart = tEnd
                qStart = qEnd
        cigar = ""
        for operation, length in operations_lengths:
            cigar += str(length) + "MIDNSHP=X"[operation]
        try:
            mapq = alignment.mapq
        except AttributeError:
            mapq = 255  # not available
        rNext = "*"
        pNext = 0
        tLen = 0
        fields = [
            qName,
            str(flag),
            rName,
            str(pos + 1),  # 1-based coordinates
            str(mapq),
            cigar,
            rNext,
            str(pNext),
            str(tLen),
            query,
            qual,
        ]
        if md is None:
            md = self.md
        if md is True:
            if query == "*":
                raise ValueError("requested MD tag with undefined sequence")
            number = 0
            ts = pos
            qs = 0
            md = ""
            for operation, length in operations_lengths:
                if operation == 0:  # M; alignment match
                    te = ts + length
                    qe = qs + length
                    for tc, qc in zip(target[ts:te], query[qs:qe]):
                        if tc == qc:
                            number += 1
                        else:
                            md += str(number) + tc
                            number = 0
                    ts = te
                    qs = qe
                elif operation == 1:  # I; insertion to the reference
                    qe = qs + length
                    qs = qe
                elif operation == 2:  # D; deletion from the reference
                    te = ts + length
                    if number:
                        md += str(number)
                        number = 0
                    md += "^" + target[ts:te]
                    ts = te
                elif operation == 3:  # N; skipped region from the reference
                    te = ts + length
                    ts = te
                elif operation == 4:  # S; soft clipping
                    qe = qs + length
                    qs = qe
                elif operation == 5:  # H; hard clipping
                    pass
                elif operation == 6:  # P; padding
                    raise NotImplementedError("padding operator is not yet implemented")
                elif operation == 7:  # =; sequence match
                    te = ts + length
                    qe = qs + length
                    number += length
                    ts = te
                    qs = qe
                elif operation == 8:  # X; sequence mismatch
                    te = ts + length
                    qe = qs + length
                    for tc in target[ts:te]:
                        md += str(number)
                        number = 0
                    ts = te
                    qs = qe
            if number:
                md += str(number)
            field = "MD:Z:%s" % md
            fields.append(field)
        try:
            score = alignment.score
        except AttributeError:
            pass
        else:
            field = "AS:i:%d" % int(round(score))
            fields.append(field)
        try:
            annotations = alignment.annotations
        except AttributeError:
            pass
        else:
            for key, value in annotations.items():
                if isinstance(value, int):
                    datatype = "i"
                    value = str(value)
                elif isinstance(value, float):
                    datatype = "f"
                    value = str(value)
                elif isinstance(value, str):
                    if len(value) == 1:
                        datatype = "A"
                    else:
                        datatype = "Z"
                elif isinstance(value, bytes):
                    datatype = "H"
                    value = "".join(map(str, value))
                elif isinstance(value, numpy.array):
                    datatype = "B"
                    if numpy.issubdtype(value.dtype, numpy.integer):
                        letter = "i"
                    elif numpy.issubdtype(value.dtype, float):
                        letter = "f"
                    else:
                        raise ValueError(
                            f"Array of incompatible data type {value.dtype} in annotation '{key}'"
                        )
                    value = ",".join(map(str, value))
                field = f"{key}:{datatype}:{value}"
                fields.append(field)
        line = "\t".join(fields) + "\n"
        return line


class AlignmentIterator(interfaces.AlignmentIterator):
    """Alignment iterator for Sequence Alignment/Map (SAM) files.

    Each line in the file contains one genomic alignment, which are loaded
    and returned incrementally.  The following columns are stored as attributes
    of the alignment:

      - flag: The FLAG combination of bitwise flags;
      - mapq: Mapping Quality (only stored if available)
      - rnext: Reference sequence name of the primary alignment of the next read
               in the alignment (only stored if available)
      - pnext: Zero-based position of the primary alignment of the next read in
               the template (only stored if available)
      - tlen: signed observed template length (only stored if available)

    Other information associated with the alignment by its tags are stored in
    the annotations attribute of each alignment.

    The sequence quality, if available, is stored as 'phred_quality' in the
    letter_annotations dictionary attribute of the query sequence record.
    """

    def __init__(self, source):
        """Create an AlignmentIterator object.

        Arguments:
         - source   - input data or file name

        """
        super().__init__(source, mode="t", fmt="SAM")
        stream = self.stream
        self.metadata = {}
        self.targets = {}
        for line in stream:
            if not line.startswith("@"):
                self.line = line
                break
            fields = line[1:].strip().split("\t")
            tag = fields[0]
            values = {}
            if tag == "SQ":
                annotations = {}
                description = None
                for field in fields[1:]:
                    key, value = field.split(":", 1)
                    assert len(key) == 2
                    if key == "SN":
                        rname = value
                    elif key == "LN":
                        length = int(value)
                    elif key == "AH":
                        annotations["alternate_locus"] = value
                    elif key == "AN":
                        annotations["names"] = value.split(",")
                    elif key == "AS":
                        annotations["assembly"] = value
                    elif key == "DS":
                        description = value
                    elif key == "M5":
                        annotations["MD5"] = value
                    elif key == "SP":
                        annotations["species"] = value
                    elif key == "TP":
                        assert value in ("linear", "circular")
                        annotations["topology"] = value
                    elif key == "UR":
                        annotations["URI"] = value
                    else:
                        annotations[key] = value
                assert rname not in self.targets
                sequence = Seq(None, length=length)
                record = SeqRecord(sequence, id=rname, annotations=annotations)
                if description is not None:
                    record.description = description
                self.targets[rname] = record
            else:
                for field in fields[1:]:
                    key, value = field.split(":", 1)
                    assert len(key) == 2
                    values[key] = value
                if tag == "HD":
                    self.metadata[tag] = values
                else:
                    if tag not in self.metadata:
                        self.metadata[tag] = []
                    self.metadata[tag].append(values)

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
            fields = line.split()
            if len(fields) < 11:
                raise ValueError(
                    "line has %d columns; expected at least 11" % len(fields)
                )
            qname = fields[0]
            flag = int(fields[1])
            rname = fields[2]
            pos = int(fields[3]) - 1
            mapq = int(fields[4])
            cigar = []
            number = ""
            for letter in fields[5]:
                if letter == "M":  # alignment match
                    operation = 0
                elif letter == "I":  # insertion to the reference
                    operation = 1
                elif letter == "D":  # deletion from the reference
                    operation = 2
                elif letter == "N":  # skipped region from the reference
                    operation = 3
                elif letter == "S":  # soft clipping
                    operation = 4
                elif letter == "H":  # hard clipping
                    operation = 5
                elif letter == "P":  # padding
                    operation = 6
                elif letter == "=":  # sequence match
                    operation = 7
                elif letter == "X":  # sequence mismatch
                    operation = 8
                else:
                    number += letter
                    continue
                cigar.append((operation, int(number)))
                number = ""
            rnext = fields[6]
            pnext = int(fields[7]) - 1
            tlen = int(fields[8])
            query = fields[9]
            qual = fields[10]
            md = None
            score = None
            annotations = {}
            column_annotations = {}
            for field in fields[11:]:
                tag, datatype, value = field.split(":", 2)
                if tag == "AS":
                    assert datatype == "i"
                    score = int(value)
                elif tag == "MD":
                    assert datatype == "Z"
                    md = value
                else:
                    if datatype == "i":
                        value = int(value)
                    elif datatype == "f":
                        value = float(value)
                    elif datatype in ("A", "Z"):  # string
                        pass
                    elif datatype == "H":
                        n = len(value)
                        value = bytes(int(value[i : i + 2]) for i in range(0, n, 2))
                    elif datatype == "B":
                        letter = value[0]
                        value = value[1:].split(",")
                        if letter in "cCsSiI":
                            dtype = int
                        elif letter == "f":
                            dtype = float
                        else:
                            raise ValueError(
                                f"Unknown number type '{letter}' in tag '{field}'"
                            )
                        value = numpy.array(value, dtype)
                    annotations[tag] = value
            if flag & 0x10:
                strand = "-"
            else:
                strand = "+"
            if pos >= 0:
                coordinates = Coordinates(cigar, pos, strand)
            else:
                coordinates = None
            if md is None:
                if rname == "*":  # unmapped
                    target = None
                else:
                    target = self.targets.get(rname)
                    if target is None:
                        if self.targets:
                            raise ValueError(
                                f"Found target {rname} missing from header"
                            )
                        target = SeqRecord(None, id=rname)
            else:
                seq = query
                target = ""
                starts = [pos]
                size = 0
                sizes = []
                number = ""
                for operation, length in cigar:
                    if operation in (0, 7, 8):  # M=X
                        pos += length
                        target += seq[:length]
                        seq = seq[length:]
                        size += length
                    elif operation in (1, 4):  # IS
                        seq = seq[length:]
                    elif operation == 2:  # D
                        pos += length
                        size += length
                        starts.append(pos)
                        sizes.append(size)
                        size = 0
                    elif operation == 3:  # N
                        pos += length
                        starts.append(pos)
                        sizes.append(size)
                        size = 0
                    elif operation == 5:  # H
                        # hard clipping (clipped sequences not present in sequence)
                        continue
                    elif operation == 6:  # P
                        raise NotImplementedError(
                            "padding operator is not yet implemented"
                        )
                sizes.append(size)
                seq = target
                target = ""
                letters = iter(md)
                number = ""
                for letter in letters:
                    if letter in "ACGTNacgtn":
                        if number:
                            number = int(number)
                            target += seq[:number]
                            seq = seq[number:]
                            number = ""
                        target += letter
                        seq = seq[1:]
                    elif letter == "^":
                        if number:
                            number = int(number)
                            target += seq[:number]
                            seq = seq[number:]
                            number = ""
                        for letter in letters:
                            if letter not in "ACGTNacgtn":
                                break
                            target += letter
                        else:
                            break
                        number = letter
                    else:
                        number += letter
                if number:
                    number = int(number)
                    target += seq[:number]
                seq = target
                target = copy.deepcopy(self.targets[rname])
                length = len(target.seq)
                data = {}
                index = 0
                for start, size in zip(starts, sizes):
                    data[start] = seq[index : index + size]
                    index += size
                target.seq = Seq(data, length=length)
            if query == "*":
                length = abs(coordinates[1, -1] - coordinates[1, 0])
                sequence = Seq(None, length=length)
            else:
                sequence = Seq(query)
                if strand == "-":
                    sequence = sequence.reverse_complement()
            query = SeqRecord(sequence, id=qname)
            if qual != "*":
                query.letter_annotations["phred_quality"] = qual
            records = [target, query]
            alignment = Alignment(records, coordinates)
            alignment.flag = flag
            if mapq != 255:
                alignment.mapq = mapq
            if rnext == "=":
                alignment.rnext = rname
            elif rnext != "*":
                alignment.rnext = rnext
            if pnext >= 0:
                alignment.pnext = pnext
            if tlen != 0:
                alignment.tlen = tlen
            if score is not None:
                alignment.score = score
            if annotations:
                alignment.annotations = annotations
            yield alignment
