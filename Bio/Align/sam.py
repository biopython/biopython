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
import numpy


from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq, reverse_complement, UndefinedSequenceError
from Bio.SeqRecord import SeqRecord


class AlignmentWriter(interfaces.AlignmentWriter):
    """Alignment file writer for the Sequence Alignment/Map (SAM) ile format."""

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
                    fields.append("DS:%s" % value)
            line = "\t".join(fields)+ "\n"
            self.stream.write(line)
        for tag, rows in metadata.items():
            if tag == "HD":  # already written
                continue
            for row in rows:
                fields = ["@" + tag]
                for key, value in row.items():
                    fields.append("%s:%s" % (key, value))
                line = "\t".join(fields)+ "\n"
                self.stream.write(line)

    def format_alignment(self, alignment):
        """Return a string with a single alignment formatted as one SAM line."""
        if not isinstance(alignment, Alignment):
            raise TypeError("Expected an Alignment object")
        target, query = alignment.sequences
        try:
            qName = query.id
        except AttributeError:
            qName = "query"
        else:
            query = query.seq
        try:
            rName = target.id
        except AttributeError:
            rName = "target"
        else:
            target = target.seq
        n1 = len(target)
        n2 = len(query)
        pos = None
        tSize = n1
        cigar = []
        coordinates = alignment.coordinates
        if coordinates[1, 0] < coordinates[1, -1]:  # mapped to forward strand
            flag = 0
            seq = query
        else:  # mapped to reverse strand
            flag = 16
            seq = reverse_complement(query, inplace=False)
            coordinates = coordinates.copy()
            coordinates[1, :] = n2 - coordinates[1, :]
        try:
            seq = bytes(seq)
        except TypeError:  # string
            pass
        except UndefinedSequenceError:
            seq = "*"
        else:
            seq = str(seq, "ASCII")
        tStart, qStart = coordinates[:, 0]
        for tEnd, qEnd in coordinates[:, 1:].transpose():
            tCount = tEnd - tStart
            qCount = qEnd - qStart
            if tCount == 0:
                length = qCount
                if pos is None or tEnd == tSize:
                    operation = "S"
                else:
                    operation = "I"
                qStart = qEnd
            elif qCount == 0:
                if tStart > 0 and tEnd < tSize:
                    length = tCount
                    operation = "D"
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
                operation = "M"
                length = tCount
            if operation is not None:
                cigar.append(str(length) + operation)
        mapQ = 255  # not available
        rNext = "*"
        pNext = 0
        tLen = 0
        try:
            qual = alignment.column_annotations["qual"]
        except (AttributeError, KeyError):
            qual = "*"
        cigar = "".join(cigar)
        fields = [
            qName,
            str(flag),
            rName,
            str(pos + 1),  # 1-based coordinates
            str(mapQ),
            cigar,
            rNext,
            str(pNext),
            str(tLen),
            seq,
            qual,
        ]
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
                        raise ValueError(f"Array of incompatible data type {value.dtype} in annotation '{key}'")
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
      - qual: ASCII of base quality plus 33 (only stored if available)

    Other information associated with the alignment by its tags are stored in
    the annotations attribute of each alignment.
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
                raise ValueError("line has %d columns; expected at least 11" % len(fields))
            qname = fields[0]
            flag = int(fields[1])
            rname = fields[2]
            pos = int(fields[3]) - 1
            mapq = int(fields[4])
            cigar = fields[5]
            rnext = fields[6]
            pnext = int(fields[7]) - 1
            tlen = int(fields[8])
            seq = fields[9]
            qual = fields[10]
            number = ""
            query_pos = 0
            coordinates = [[pos, query_pos]]
            for letter in cigar:
                if letter in "M=X":
                    number = int(number)
                    pos += number
                    query_pos += number
                    coordinates.append([pos, query_pos])
                    number = ""
                elif letter in "IS":
                    number = int(number)
                    query_pos += number
                    coordinates.append([pos, query_pos])
                    number = ""
                elif letter in "DN":
                    number = int(number)
                    pos += number
                    coordinates.append([pos, query_pos])
                    number = ""
                elif letter == "H":
                    # hard clipping (clipped sequences not present in sequence)
                    number = int(number)
                    number = ""
                elif letter == "P":
                    number = int(number)
                    raise NotImplementedError("padding operator is not yet implemented")
                    number = ""
                else:
                    number += letter
            coordinates = numpy.array(coordinates).transpose()
            target = self.targets[rname]
            if seq == "*":
                sequence = Seq(None, length=query_pos)
            else:
                sequence = Seq(seq)
                if flag & 0x10:
                    sequence = sequence.reverse_complement()
            if flag & 0x10:
                coordinates[1, :] = query_pos - coordinates[1, :]
            query = SeqRecord(sequence, id=qname)
            records = [target, query]
            alignment = Alignment(records, coordinates)
            annotations = {}
            column_annotations = {}
            alignment.flag = flag
            if mapq != 255:
                alignment.mapq = mapq
            if rnext == "=":
                alignment.rnext = rname
            elif rnext != "*":
                alignment.rnext = rnext
            if pnext >= 0:
                alignment.pnext = pnext
            if tlen > 0:
                alignment.tlen = tlen
            if qual != "*":
                column_annotations["quality"] = qual
            for field in fields[11:]:
                tag, datatype, value = field.split(":", 2)
                if tag == "AS":
                    assert datatype == "i"
                    alignment.score = int(value)
                else:
                    if datatype == "i":
                        value = int(value)
                    elif datatype == "f":
                        value = float(value)
                    elif datatype in ("A", "Z"):  # string
                        pass
                    elif datatype == "H":
                        n = len(value)
                        value = bytes(int(value[i: i + 2]) for i in range(0, n, 2))
                    elif datatype == "B":
                        letter = value[0]
                        value = value[1:].split(",")
                        if letter in "cCsSiI":
                            dtype = int
                        elif letter == "f":
                            dtype = float
                        else:
                            raise ValueError(f"Unknown number type '{letter}' in tag '{field}'")
                        value = numpy.array(value, dtype)
                    annotations[tag] = value
            if annotations:
                alignment.annotations = annotations
            if column_annotations:
                alignment.column_annotations = column_annotations
            yield alignment
